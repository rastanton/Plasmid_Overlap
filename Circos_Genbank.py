from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch
import glob

def Dash_Replace(input_string):
    """Replaces - w/ _"""
    Out = ''
    for entry in input_string:
        if entry == '-':
            Out = Out + '_'
        else:
            Out = Out + entry
    return Out

def Dash_Replace_List(input_list):
    """Repaces - with _ for all names in a list"""
    Out = []
    for entry in input_list:
        Name = Dash_Replace(entry)
        Out.append(entry)
    return Out

def Repeat_Remover(any_list):
    """Removes repeats for any list"""
    new_list = []
    for items in any_list:
        if (items in new_list) == False:
            new_list.append(items)
    return new_list

def GenBank_Feature_Names(input_genbank):
    f = open(input_genbank, 'r')
    Features = []
    Add = 0
    for line in f:
        Info = line.split()
        if Info[0] == 'ORIGIN':
            Features = Repeat_Remover(Features)
            return Features
        elif Info[0] == 'FEATURES':
            Add = 1
        elif Add == 1:
            Features.append(Info[0])

def GenBank_Figure(input_gb, match_name_list, match_color_list, title, output_file):
    """Makes a circos GB match figure"""
    gbk = Genbank(input_gb)
    circos = Circos(sectors={gbk.name: gbk.range_size})
    circos.text(title, size=12, r=20)
    sector = circos.get_sector(gbk.name)
    major_ticks_interval = 10000
    minor_ticks_interval = 1000
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 3:.1f} kb"
        )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)
    f_cds_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="slategrey")
    r_cds_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="slategrey")
    AMR_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    AMR_f_track.genomic_features(gbk.extract_features("AMR", target_strand=1), plotstyle="arrow", fc="red")
    AMR_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    AMR_r_track.genomic_features(gbk.extract_features("AMR", target_strand=-1), plotstyle="arrow", fc="red")
    PF_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    PF_f_track.genomic_features(gbk.extract_features("PF", target_strand=1), plotstyle="arrow", fc="blue")
    PF_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    PF_r_track.genomic_features(gbk.extract_features("PF", target_strand=-1), plotstyle="arrow", fc="blue")
    Names = Dash_Replace_List(match_name_list)
    Position = 92
    for entry in range(len(Names)):
        Track = "ID_" + Names[entry] + "_track"
        Track = sector.add_track(((Position - 5), Position), r_pad_ratio=0.1)
        Position = Position - 5
        Track.genomic_features(gbk.extract_features(match_name_list[entry]), fc=match_color_list[entry])
    handles = []
    for entry in range(len(match_name_list)):
        handles.append(Patch(color=match_color_list[entry], label=match_name_list[entry]))
    fig = circos.plotfig()
    _ = fig.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=12)
    fig.savefig(output_file)

def GenBank_Figure_No_Label(input_gb, match_name_list, match_color_list, output_file):
    """Makes a circos GB match figure"""
    gbk = Genbank(input_gb)
    circos = Circos(sectors={gbk.name: gbk.range_size})
    #circos.text(title, size=12, r=20)
    sector = circos.get_sector(gbk.name)
    major_ticks_interval = 10000
    minor_ticks_interval = 1000
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 3:.1f} kb"
        )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)
    f_cds_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="slategrey")
    r_cds_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="slategrey")
    AMR_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    AMR_f_track.genomic_features(gbk.extract_features("AMR", target_strand=1), plotstyle="arrow", fc="red")
    AMR_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    AMR_r_track.genomic_features(gbk.extract_features("AMR", target_strand=-1), plotstyle="arrow", fc="red")
    PF_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    PF_f_track.genomic_features(gbk.extract_features("PF", target_strand=1), plotstyle="arrow", fc="blue")
    PF_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    PF_r_track.genomic_features(gbk.extract_features("PF", target_strand=-1), plotstyle="arrow", fc="blue")
    Names = Dash_Replace_List(match_name_list)
    Names = Dash_Replace_List(match_name_list)
    Position = 92
    for entry in range(len(Names)):
        Track = "ID_" + Names[entry] + "_track"
        Track = sector.add_track(((Position - 5), Position), r_pad_ratio=0.1)
        Position = Position - 5
        Track.genomic_features(gbk.extract_features(match_name_list[entry]), fc=match_color_list[entry])
    #handles = []
    #for entry in range(len(match_name_list)):
    #    handles.append(Patch(color=match_color_list[entry], label=match_name_list[entry]))
    fig = circos.plotfig()
    #_ = fig.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=12)
    fig.savefig(output_file)

def GenBank_Figure_No_Label_Test(input_gb, match_name_list, match_color_list, title, output_file):
    """Makes a circos GB match figure"""
    gbk = Genbank(input_gb)
    circos = Circos(sectors={gbk.name: gbk.range_size})
    #circos.text(title, size=12, r=20)
    sector = circos.get_sector(gbk.name)
    major_ticks_interval = 10000
    minor_ticks_interval = 1000
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 3:.1f} kb"
        )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)
    f_cds_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="slategrey")
    r_cds_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="slategrey")
    AMR_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    AMR_f_track.genomic_features(gbk.extract_features("AMR", target_strand=1), plotstyle="arrow", fc="red")
    AMR_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    AMR_r_track.genomic_features(gbk.extract_features("AMR", target_strand=-1), plotstyle="arrow", fc="red")
    PF_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    PF_f_track.genomic_features(gbk.extract_features("PF", target_strand=1), plotstyle="arrow", fc="blue")
    PF_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    PF_r_track.genomic_features(gbk.extract_features("PF", target_strand=-1), plotstyle="arrow", fc="blue")
    Names = Dash_Replace_List(match_name_list)
    Position = 92
##    for entry in range(len(Names)):
##        Track = "ID_" + Names[entry] + "_track"
##        Track = sector.add_track(((Position - 5), Position), r_pad_ratio=0.1)
##        Position = Position - 5
##        Track.genomic_features(gbk.extract_features(match_name_list[entry]), fc=match_color_list[entry])
    #handles = []
    #for entry in range(len(match_name_list)):
    #    handles.append(Patch(color=match_color_list[entry], label=match_name_list[entry]))
    fig = circos.plotfig()
    #_ = fig.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=12)
    fig.savefig(output_file)

def GenBank_Figure_No_Label_2(input_gb, match_name_list, match_color_list, output_file):
    """Makes a circos GB match figure"""
    gbk = Genbank(input_gb)
    circos = Circos(sectors={gbk.name: gbk.range_size})
    #circos.text(title, size=12, r=20)
    sector = circos.get_sector(gbk.name)
    major_ticks_interval = 10000
    minor_ticks_interval = 1000
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 3:.1f} kb"
        )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)
##    f_cds_track = sector.add_track((95, 97), r_pad_ratio=0.1)
##    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="slategrey")
##    r_cds_track = sector.add_track((93, 95), r_pad_ratio=0.1)
##    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="slategrey")
    AMR_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    AMR_f_track.genomic_features(gbk.extract_features("AMR", target_strand=1), plotstyle="arrow", fc="red")
    AMR_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    AMR_r_track.genomic_features(gbk.extract_features("AMR", target_strand=-1), plotstyle="arrow", fc="red")
    PF_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    PF_f_track.genomic_features(gbk.extract_features("PF", target_strand=1), plotstyle="arrow", fc="blue")
    PF_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    PF_r_track.genomic_features(gbk.extract_features("PF", target_strand=-1), plotstyle="arrow", fc="blue")
    Names = Dash_Replace_List(match_name_list)
    Position = 92
    for entry in range(len(Names)):
        Track = "ID_" + Names[entry] + "_track"
        Track = sector.add_track(((Position - 5), Position), r_pad_ratio=0.1)
        Position = Position - 5
        Track.genomic_features(gbk.extract_features(match_name_list[entry]), fc=match_color_list[entry])
    #handles = []
    #for entry in range(len(match_name_list)):
    #    handles.append(Patch(color=match_color_list[entry], label=match_name_list[entry]))
    fig = circos.plotfig()
    #_ = fig.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=12)
    fig.savefig(output_file)

def GenBank_Figure_No_Label_Small(input_gb, match_name_list, match_color_list, output_file):
    """Makes a circos GB match figure"""
    gbk = Genbank(input_gb)
    circos = Circos(sectors={gbk.name: gbk.range_size})
    #circos.text(title, size=12, r=20)
    sector = circos.get_sector(gbk.name)
    major_ticks_interval = 10000
    minor_ticks_interval = 1000
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{v/ 10 ** 3:.1f} kb"
        )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)
##    f_cds_track = sector.add_track((95, 97), r_pad_ratio=0.1)
##    f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), plotstyle="arrow", fc="slategrey")
##    r_cds_track = sector.add_track((93, 95), r_pad_ratio=0.1)
##    r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), plotstyle="arrow", fc="slategrey")
    AMR_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    AMR_f_track.genomic_features(gbk.extract_features("AMR", target_strand=1), plotstyle="arrow", fc="red")
    AMR_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    AMR_r_track.genomic_features(gbk.extract_features("AMR", target_strand=-1), plotstyle="arrow", fc="red")
    PF_f_track = sector.add_track((95, 97), r_pad_ratio=0.1)
    PF_f_track.genomic_features(gbk.extract_features("PF", target_strand=1), plotstyle="arrow", fc="blue")
    PF_r_track = sector.add_track((93, 95), r_pad_ratio=0.1)
    PF_r_track.genomic_features(gbk.extract_features("PF", target_strand=-1), plotstyle="arrow", fc="blue")
    Names = Dash_Replace_List(match_name_list)
    Position = 92
    for entry in range(len(Names)):
        Track = "ID_" + Names[entry] + "_track"
        Track = sector.add_track(((Position - 3), Position), r_pad_ratio=0.1)
        Position = Position - 3
        Track.genomic_features(gbk.extract_features(match_name_list[entry]), fc=match_color_list[entry])
    #handles = []
    #for entry in range(len(match_name_list)):
    #    handles.append(Patch(color=match_color_list[entry], label=match_name_list[entry]))
    fig = circos.plotfig()
    #_ = fig.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=12)
    fig.savefig(output_file)
