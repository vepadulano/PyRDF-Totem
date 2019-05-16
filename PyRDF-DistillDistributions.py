import PyRDF
import ROOT
from ROOT import TCanvas

# from pyspark import SparkContext

# context = SparkContext.getOrCreate()
# context.stop()

PyRDF.use("spark", {'npartitions':2})

PyRDF.include_headers("common_definitions.h")
PyRDF.include_headers("parameters_global.h")
PyRDF.include_headers("common_algorithms.h")
PyRDF.include_headers("parameters.h")
PyRDF.include_headers("common.h")
PyRDF.include_headers("initialize.h")

#####################################################
# Prepare input data
#####################################################

# Branches clasified by diagonal
diagonals = {
    # return a tuple: ([left] verticals in 45, [right] verticals in 56))
    "d45b_56t"  : (["track_rp_5", "track_rp_21", "track_rp_25"],
                   ["track_rp_104", "track_rp_120", "track_rp_124"]),
    "ad45b_56b" : (["track_rp_5", "track_rp_21", "track_rp_25"],
                   ["track_rp_105", "track_rp_121", "track_rp_125"]),
    "d45t_56b"  : (["track_rp_4", "track_rp_20", "track_rp_24"],
                   ["track_rp_105", "track_rp_121", "track_rp_125"]),
    "ad45t_56t" : (["track_rp_4", "track_rp_20", "track_rp_24"],
                   ["track_rp_104", "track_rp_120", "track_rp_124"])
}

# Data source files
DS = {
    'DS0'   : '4495',
    'DS1'   : '4495',
    'DS2'   : '4496',
    'DS3'   : '4499',
    'DS4'   : '4505',
    'DS5'   : '4509',
    'DS6'   : '4510',
    'DS7'   : '4511'
}

# Select branches
selected_diagonal = "d45b_56t"
selected_ds = 'DS0'

rp_left, rp_right = diagonals[selected_diagonal]

# Columns per branch
attributes = ['valid', 'x', 'y']

full_branches = ["{}.{}".format(c,a) for a in attributes for c in rp_left+rp_right ]

# Split left and right branch on valid, x and y
valids = [ "(unsigned int) {}".format(v) for v in full_branches[0:6]]
xs     = full_branches[6:12]
ys     = full_branches[12:18]

print("Selected branches: \n" + "\n\t".join(full_branches))

# Filter to select valid branches
filter_code = """({0}.valid + {1}.valid + {2}.valid ) >= 2 &&
({3}.valid + {4}.valid + {5}.valid ) >= 2
""".format(*(rp_left+rp_right))

# Original analysis skips some input files, hence we use input_files*.txt files
# to specify the subset of processed input files.
source_file       = "input_files_{}.txt".format(selected_ds)

# Name of the ROOT Tree
input_ntuple_name = "TotemNtuple"

# Remote path in EOS
#prefix            = "root://eostotem//eos/totem/data/cmstotem/2015/90m/Totem/Ntuple/version2/{}/".format(DS[selected_ds])

# Read from local disk
prefix            = "/data/TotemData/{}/".format(DS[selected_ds])

# Construct list of prefix + input files
input_files       = [prefix + line.rstrip('\r\n') for line in open(source_file)]

# Create chain of files
c = ROOT.TChain("TotemNtuple")
for filename in input_files:
  c.Add(filename)


############################
# Analysis initialization
############################

import ROOT

# Load C++ headers
PyRDF.include_headers("common_definitions.h")
PyRDF.include_headers("parameters_global.h")
PyRDF.include_headers("common_algorithms.h")
PyRDF.include_headers("parameters.h")
PyRDF.include_headers("common.h")
PyRDF.include_headers("initialize.h")


detailsLevel         = 0 # 0: no details, 1: some details, >= 2 all details
overrideCutSelection = False # whether the default cut selection should be overriden by the command-line selection
cutSelectionString   = None
outputDir            = "."
inputDir             = "."
input_n_si           = 4.0
time_group_divisor   = 0
time_group_remainder = 0
event_group_divisor  = 0
event_group_index    = 0
evIdxStep            = 1
maxTaggedEvents      = 0

M_PI = 3.14159265358979323846264338328      # Pi


def myinit():
    import ROOT
    ROOT.initialize()

PyRDF.initialize(myinit)

# TODO
#f = getattr(ROOT, "initialize")
#f(*args)

# alignment init
# for i,_ in  enumerate(ROOT.alignmentSources):
#     ROOT.alignmentSources[i].Init()

# get time-dependent corrections
corrg_pileup = None
if ROOT.anal.use_pileup_efficiency_fits:
    path = inputDir + "/pileup_fit_combined.root"
    puF = ROOT.TFile.Open(path)
    if not os.path.exists(puF):
        print("ERROR: pile-up correction file `%s' cannot be opened.\n" % path);
    if diagonal == "d45b_56t":
        #corrg_pileup = (TGraph *) puF.Get("45b_56t/dgn");
        corrg_pileup = puF.Get("45b_56t/dgn")
    if diagonal == "d45t_56b":
        #corrg_pileup = (TGraph *) puF.Get("45b_56t/dgn");
        corrg_pileup = puF.Get("45t_56b/dgn")

# get th_y* dependent efficiency correction
f_3outof4_efficiency_L_F = None
f_3outof4_efficiency_L_N = None
f_3outof4_efficiency_R_N = None
f_3outof4_efficiency_R_F = None

if ROOT.anal.use_3outof4_efficiency_fits:
    path = inputDir + "/eff3outof4_details_fit_old.root"
    effFile = ROOT.TFile.Open(path)
    if (os.path.exists(effFile)):
        print("ERROR: 3-out-of-4 efficiency file `%s' cannot be opened.\n" % path)

    diagonal = selected_diagonal;
    f_3outof4_efficiency_L_F = effFile.Get(diagonal + "/L_F/fit")
    f_3outof4_efficiency_L_N = effFile.Get(diagonal + "/L_N/fit")
    f_3outof4_efficiency_R_N = effFile.Get(diagonal + "/R_N/fit")
    f_3outof4_efficiency_R_F = effFile.Get(diagonal + "/R_F/fit")

    print("\n>> using 3-out-of-4 fits: %s, %s, %s, %s\n" % (f_3outof4_efficiency_L_F,
                                                            f_3outof4_efficiency_L_N,
                                                            f_3outof4_efficiency_R_N,
                                                            f_3outof4_efficiency_R_F))

# book metadata histograms
# ROOT.timestamp_bins = ROOT.timestamp_max - ROOT.timestamp_min + 1.

# Create bh_t_* hists
bh_t_Nev_before = dict()
bh_t_Nev_after_no_corr = dict()
bh_t_before = dict()
bh_t_after = dict()
bh_t_after_no_corr = dict()
bp_t_phi_corr = dict()
bp_t_full_corr = dict()

# binnings
binnings = ["ub", "ob-1-10-0.2", "ob-1-30-0.2"]

binning_setup = dict()
for b in binnings:
    binning = ROOT.BuildBinningRDF(ROOT.anal, b)
    binning_setup[b] = binning

# zero counters
n_ev_full = 0;
n_ev_cut = dict();
for ci in range(ROOT.anal.N_cuts):
    n_ev_cut[ci] = 0

th_min = 1E100;
th_y_L_min = +1E100; th_y_R_min = +1E100

N_anal=0; N_anal_zeroBias=0;
N_zeroBias_el=0; N_zeroBias_el_RP_trig=0;
N_4outof4=0; N_el=0;
N_el_T2trig=0; N_4outof4_T2trig=0;
N_el_raw=0;

########################
# Initialize RDataFrame
########################

df = PyRDF.RDataFrame("TotemNtuple", "../TotemNTuple_9883.040.ntuple.root")

# Distill

# Filter and define output branches
distilled = df.Filter(filter_code)\
            .Define("v_L_1_F", valids[0])\
            .Define("v_L_2_N", valids[1])\
            .Define("v_L_2_F", valids[2])\
            .Define("v_R_1_F", valids[3])\
            .Define("v_R_2_N", valids[4])\
            .Define("v_R_2_F", valids[5])\
            .Define("x_L_1_F", xs[0])\
            .Define("x_L_2_N", xs[1])\
            .Define("x_L_2_F", xs[2])\
            .Define("x_R_1_F", xs[3])\
            .Define("x_R_2_N", xs[4])\
            .Define("x_R_2_F", xs[5])\
            .Define("y_L_1_F", ys[0])\
            .Define("y_L_2_N", ys[1])\
            .Define("y_L_2_F", ys[2])\
            .Define("y_R_1_F", ys[3])\
            .Define("y_R_2_N", ys[4])\
            .Define("y_R_2_F", ys[5])\
            .Define("timestamp",    "(unsigned int) (event_info.timestamp - 1444860000)")\
            .Define("run_num",      "(unsigned int) event_info.run_no")\
            .Define("bunch_num",    "trigger_data.bunch_num")\
            .Define("event_num",    "trigger_data.event_num")\
            .Define("trigger_num",  "trigger_data.trigger_num")\
            .Define("trigger_bits", "trigger_data.input_status_bits")

# Distributions

#########################################################
###### FILTER, BUILD HISTOGRAMS - START EVENT LOOP ######
#########################################################

rdf = distilled

# Initial cuts
f1 = rdf.Filter("! SkipTime( timestamp )", 'check time - selected')

if time_group_divisor != 0:
    f1 = f1.Filter("! SkipTimeInterval( timestamp, %s, %s )".format(time_group_divisor, time_group_remainder),
                    'time interval')

# Diagonal cut (L831)
f2 = f1.Filter("v_L_2_F && v_L_2_N && v_R_2_F && v_R_2_N", 'allDiagonalRPs')

model = ("h_timestamp_dgn", ";timestamp;rate   (Hz)",
         int(ROOT.timestamp_bins), ROOT.timestamp_min-0.5, ROOT.timestamp_max+0.5)
h_timestamp_dgn = f2.Histo1D(model, "timestamp")

# Not cut for this filter in original code
f_zerobias = f2.Filter("! ((trigger_bits & 512) != 0)", 'zero_bias_event')

xs = ["x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F"]
ys = ["y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F"]

# Apply fine alignment (L 852)
r2 = f2.Define("h_al", "ApplyFineAlignment( timestamp ,{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(*(xs+ys) ))         .Define("h_al_x_L_1_F", "h_al.L_1_F.x")         .Define("h_al_x_L_2_N", "h_al.L_2_N.x")         .Define("h_al_x_L_2_F", "h_al.L_2_F.x")         .Define("h_al_y_L_1_F", "h_al.L_1_F.y")         .Define("h_al_y_L_2_N", "h_al.L_2_N.y")         .Define("h_al_y_L_2_F", "h_al.L_2_F.y")         .Define("h_al_x_R_1_F", "h_al.R_1_F.x")         .Define("h_al_x_R_2_N", "h_al.R_2_N.x")         .Define("h_al_x_R_2_F", "h_al.R_2_F.x")         .Define("h_al_y_R_1_F", "h_al.R_1_F.y")         .Define("h_al_y_R_2_N", "h_al.R_2_N.y")         .Define("h_al_y_R_2_F", "h_al.R_2_F.y")
# fill pre-selection histograms
al_nosel_models = [
    ("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.),
    ("h_y_L_2_N_vs_x_L_2_N_al_nosel", ";x^{L,2,N};y^{L,2,N}", 150, -15., 15., 300, -30., +30.),
    ("h_y_L_2_F_vs_x_L_2_F_al_nosel", ";x^{L,2,F};y^{L,2,F}", 150, -15., 15., 300, -30., +30.),
    ("h_y_R_1_F_vs_x_R_1_F_al_nosel", ";x^{R,1,F};y^{R,1,F}", 150, -15., 15., 300, -30., +30.),
    ("h_y_R_2_N_vs_x_R_2_N_al_nosel", ";x^{R,2,N};y^{R,2,N}", 150, -15., 15., 300, -30., +30.),
    ("h_y_R_2_F_vs_x_R_2_F_al_nosel", ";x^{R,2,F};y^{R,2,F}", 150, -15., 15., 300, -30., +30.)
]

h_y_L_1_F_vs_x_L_1_F_al_nosel = r2.Histo2D(al_nosel_models[0], "h_al_x_L_1_F", "h_al_y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_al_nosel = r2.Histo2D(al_nosel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_al_nosel = r2.Histo2D(al_nosel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_al_nosel = r2.Histo2D(al_nosel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_al_nosel = r2.Histo2D(al_nosel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_al_nosel = r2.Histo2D(al_nosel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F")

# run reconstruction
### kinematics struct
r3 = r2.Define("kinematics", 'DoReconstruction( h_al )')

r4 =  r3.Define("k_th_x_R",        "kinematics.th_x_R")               .Define("k_th_y_R",        "kinematics.th_y_R")               .Define("k_th_x_L",        "kinematics.th_x_L")               .Define("k_th_y_L",        "kinematics.th_y_L")               .Define("k_th_x",          "kinematics.th_x")                 .Define("k_th_y",          "kinematics.th_y")                 .Define("minus_k_th_y",    "- kinematics.th_y")               .Define("k_vtx_x",         "kinematics.vtx_x")                .Define("k_vtx_x_L",       "kinematics.vtx_x_L")              .Define("k_vtx_x_R",       "kinematics.vtx_x_R")              .Define("k_vtx_y",         "kinematics.vtx_y")                .Define("k_vtx_y_L",       "kinematics.vtx_y_L")              .Define("k_vtx_y_R",       "kinematics.vtx_y_R")              .Define("k_th_y_L_F",      "kinematics.th_y_L_F")             .Define("k_th_y_L_N",      "kinematics.th_y_L_N")             .Define("k_th_y_R_F",      "kinematics.th_y_R_F")             .Define("k_th_y_R_N",      "kinematics.th_y_R_N")             .Define("k_th_x_diffLR",   "kinematics.th_x_R - kinematics.th_x_L")                 .Define("k_th_y_diffLR",   "kinematics.th_y_R - kinematics.th_y_L")                 .Define("k_th_x_diffLF",   "kinematics.th_x_L - kinematics.th_x")                   .Define("k_th_x_diffRF",   "kinematics.th_x_R - kinematics.th_x")                   .Define("k_th_y_L_diffNF", "kinematics.th_y_L_F - kinematics.th_y_L_N")             .Define("k_th_y_R_diffNF", "kinematics.th_y_R_F - kinematics.th_y_R_N")             .Define("k_vtx_x_diffLR",  "kinematics.vtx_x_R - kinematics.vtx_x_L")               .Define("k_vtx_y_diffLR",  "kinematics.vtx_y_R - kinematics.vtx_y_L")               .Define("k_t",             "kinematics.t")                                          .Define("k_th",            "kinematics.th")                                         .Define("k_phi",           "kinematics.phi")

# cut evaluation
r5 = r4.Define("cutdata", "EvaluateCutsRDF( h_al, kinematics, anal )")

# Elastic cut
f4 = r5.Filter("cutdata.select", "elastic cut")

# Fill raw histograms
model = ("h_timestamp_sel", ";timestamp;rate   (Hz)",
         int(ROOT.timestamp_bins), ROOT.timestamp_min-0.5, ROOT.timestamp_max+0.5)
h_timestamp_sel = f4.Histo1D(model, "timestamp");

# fill histograms
noal_sel_models = list(map(ROOT.ROOT.RDF.TH2DModel, [
    ("h_y_L_1_F_vs_x_L_1_F_noal_sel", ";x^{L,1,F};y^{L,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_N_vs_x_L_2_N_noal_sel", ";x^{L,2,N};y^{L,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_F_vs_x_L_2_F_noal_sel", ";x^{L,2,F};y^{L,2,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_1_F_vs_x_R_1_F_noal_sel", ";x^{R,1,F};y^{R,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_N_vs_x_R_2_N_noal_sel", ";x^{R,2,N};y^{R,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_F_vs_x_R_2_F_noal_sel", ";x^{R,2,F};y^{R,2,F}", 100, -3., +3., 300, -30., +30.)
]))

h_y_L_1_F_vs_x_L_1_F_noal_sel = f4.Histo2D(noal_sel_models[0], "x_L_1_F", "y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_noal_sel = f4.Histo2D(noal_sel_models[1], "x_L_2_N", "y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_noal_sel = f4.Histo2D(noal_sel_models[2], "x_L_2_F", "y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_noal_sel = f4.Histo2D(noal_sel_models[3], "x_R_1_F", "y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_noal_sel = f4.Histo2D(noal_sel_models[4], "x_R_2_N", "y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_noal_sel = f4.Histo2D(noal_sel_models[5], "x_R_2_F", "y_R_2_F")

al_sel_models = list(map(ROOT.ROOT.RDF.TH2DModel, [
    ("h_y_L_1_F_vs_x_L_1_F_al_sel", ";x^{L,1,F};y^{L,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_N_vs_x_L_2_N_al_sel", ";x^{L,2,N};y^{L,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_F_vs_x_L_2_F_al_sel", ";x^{L,2,F};y^{L,2,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_1_F_vs_x_R_1_F_al_sel", ";x^{R,1,F};y^{R,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_N_vs_x_R_2_N_al_sel", ";x^{R,2,N};y^{R,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_F_vs_x_R_2_F_al_sel", ";x^{R,2,F};y^{R,2,F}", 100, -3., +3., 300, -30., +30.)
]))

h_y_L_1_F_vs_x_L_1_F_al_sel = f4.Histo2D(al_sel_models[0], "h_al_x_L_1_F", "h_al_y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_al_sel = f4.Histo2D(al_sel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_al_sel = f4.Histo2D(al_sel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_al_sel = f4.Histo2D(al_sel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_al_sel = f4.Histo2D(al_sel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_al_sel = f4.Histo2D(al_sel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F")

# Line 1157 (k.th_x_R - k.th_x_L)
#           (k.th_y_R - k.th_y_L)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("th_x_diffLR", ";#theta_{x}^{R} - #theta_{x}^{L}", 1000, -500E-6, +500E-6),
    ("th_y_diffLR", ";#theta_{y}^{R} - #theta_{y}^{L}", 500, -50E-6, +50E-6)
]))
th_x_diffLR = f4.Histo1D(models[0], "k_th_x_diffLR")
th_y_diffLR = f4.Histo1D(models[1], "k_th_y_diffLR")

# Line 1160 (k.th_x_L - k.th_x)
#           (k.th_x_R - k.th_x)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("th_x_diffLF", ";#theta_{x}^{L} - #theta_{x}", 400, -200E-6, +200E-6),
    ("th_x_diffRF", ";#theta_{x}^{R} - #theta_{x}", 400, -200E-6, +200E-6)
]))
th_x_diffLF = f4.Histo1D(models[0], "k_th_x_diffLF")
th_x_diffRF = f4.Histo1D(models[1], "k_th_x_diffRF")

# Line 1163 (k.th_x, k.th_x_R - k.th_x_L)
#           (k.th_y, k.th_y_R - k.th_y_L)
#           (k.vtx_x, k.th_x_R - k.th_x_L)
models = list(map(ROOT.ROOT.RDF.TH2DModel, [
    ("h_th_x_diffLR_vs_th_x", ";#theta_{x};#theta_{x}^{R} - #theta_{x}^{L}", 100, -300E-6, +300E-6, 120, -120E-6, +120E-6),
    ("h_th_y_diffLR_vs_th_y", ";#theta_{y};#theta_{y}^{R} - #theta_{y}^{L}", 100, -500E-6, +500E-6, 120, -120E-6, +120E-6),
    ("h_th_x_diffLR_vs_vtx_x", ";vtx_{x};#theta_{x}^{R} - #theta_{x}^{L}", 100, -300E-3, +300E-3, 120, -120E-6, +120E-6)
]))
h_th_x_diffLR_vs_th_x  = f4.Histo2D(models[0], "k_th_x", "k_th_x_diffLR")
h_th_y_diffLR_vs_th_y  = f4.Histo2D(models[1], "k_th_y", "k_th_y_diffLR")
h_th_x_diffLR_vs_vtx_x = f4.Histo2D(models[2], "k_vtx_x", "k_th_x_diffLR")


# Line 1188 (k.th_x_L, k.th_y_L)
#           (k.th_x_R, k.th_y_R)
#           (k.th_x, k.th_y)
models = list(map(ROOT.ROOT.RDF.TH2DModel, [
    ("h_th_y_L_vs_th_x_L", ";#theta_{x}^{L};#theta_{y}^{L}", 100, -115E-6, +11E-5, 100, 22E-6, +102E-6),
    ("h_th_y_R_vs_th_x_R", ";#theta_{x}^{R};#theta_{y}^{R}", 100, -125E-6, +12E-5, 100, 27E-6, +102E-6),
    ("h_th_y_vs_th_x", ";#theta_{x};#theta_{y}", 100, -300E-6, +300E-6, 100, -150E-6, +150E-6)
]))
h_th_y_L_vs_th_x_L = f4.Histo2D(models[0], "k_th_x_L", "k_th_y_L")
h_th_y_R_vs_th_x_R = f4.Histo2D(models[1], "k_th_x_R", "k_th_y_R")
h_th_y_vs_th_x     = f4.Histo2D(models[2], "k_th_x", "k_th_y")

# Line 1199 (k.th_y_R, k.th_y_L)
model = ROOT.ROOT.RDF.TH2DModel("h_th_y_L_vs_th_y_R", ";#theta_{y}^{R};#theta_{y}^{L}",
                                300, -150E-6, +150E-6, 300, -150E-6, +150E-6)
h_th_y_L_vs_th_y_R = f4.Histo2D(model, "k_th_y_R", "k_th_y_L")

# Line 1203: (k.th_x)
#            (k.th_y)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_th_x", ";#theta_{x}", 250, -500E-6, +500E-6),
    ("h_th_y", ";#theta_{y}", 250, -500E-6, +500E-6)
]))
h_th_x     = f4.Histo1D(models[0], "k_th_x")
h_th_y     = f4.Histo1D(models[1], "k_th_y")

# Line 1205: (-k.th_y)
model = ("h_th_y_flipped", ";#theta_{y}", 250, -500E-6, +500E-6)
h_th_y_flipped = f4.Histo1D(model, "minus_k_th_y")

# Line 1207: (k.th_x_L)
#            (k.th_x_R)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_th_x_L", ";#theta_{x}^{L}", 250, -500E-6, +500E-6),
    ("h_th_x_R", ";#theta_{x}^{R}", 250, -500E-6, +500E-6)
]))
h_th_x_L   = f4.Histo1D(models[0], "k_th_x_L")
h_th_x_R   = f4.Histo1D(models[1], "k_th_x_R")

# Line 1210: (k.th_y_L)
#            (k.th_y_R)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_th_y_L", ";#theta_{y}^{L}", 250, -500E-6, +500E-6),
    ("h_th_y_R", ";#theta_{y}^{R}", 250, -500E-6, +500E-6)
]))
h_th_y_L   = f4.Histo1D(models[0], "k_th_y_L")
h_th_y_R   = f4.Histo1D(models[1], "k_th_y_R")

# Line 1213: (k.th_y_L_F)
#            (k.th_y_L_N)
#            (k.th_y_R_N)
#            (k.th_y_R_F)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_th_y_L_F", ";#theta_{y}^{L_F}", 250, -500E-6, +500E-6),
    ("h_th_y_L_N", ";#theta_{y}^{L_N}", 250, -500E-6, +500E-6),
    ("h_th_y_R_N", ";#theta_{y}^{R_N}", 250, -500E-6, +500E-6),
    ("h_th_y_R_F", ";#theta_{y}^{R_F}", 250, -500E-6, +500E-6)
]))
h_th_y_L_F = f4.Histo1D(models[0], "k_th_y_L_F")
h_th_y_L_N = f4.Histo1D(models[1], "k_th_y_L_N")
h_th_y_R_N = f4.Histo1D(models[2], "k_th_y_R_N")
h_th_y_R_F = f4.Histo1D(models[3], "k_th_y_R_F")


# fill vertex histograms

# Line 1220 (k.vtx_x)
#           (k.vtx_x_L)
#           (k.vtx_x_R)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_vtx_x", ";x^{*}", 100, -0.5, +0.5),
    ("h_vtx_x_L", ";x^{*,L}", 100, -0.5, +0.5),
    ("h_vtx_x_R", ";x^{*,R}", 100, -0.5, +0.5)
]))
h_vtx_x    = f4.Histo1D(models[0], "k_vtx_x")
h_vtx_x_L  = f4.Histo1D(models[1], "k_vtx_x_L")
h_vtx_x_R  = f4.Histo1D(models[2], "k_vtx_x_R")

# Line 1224 (k.vtx_y)
#           (k.vtx_y_L)
#           (k.vtx_y_R)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_vtx_y", ";y^{*}", 100, -0.5, +0.5),
    ("h_vtx_y_L", ";y^{*,L}", 100, -0.5, +0.5),
    ("h_vtx_y_R", ";y^{*,R}", 100, -0.5, +0.5)
]))
h_vtx_y    = f4.Histo1D(models[0], "k_vtx_y")
h_vtx_y_L  = f4.Histo1D(models[1], "k_vtx_y_L")
h_vtx_y_R  = f4.Histo1D(models[2], "k_vtx_y_R")

# Line 1228:
#            (k.vtx_x_R, k.vtx_x_L)
#            (k.vtx_y_R, k.vtx_y_L)
models = list(map(ROOT.ROOT.RDF.TH2DModel, [
    ("h_vtx_x_L_vs_vtx_x_R", ";x^{*,R};x^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5),
    ("h_vtx_y_L_vs_vtx_y_R", ";y^{*,R};y^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5)
]))
h_vtx_x_L_vs_vtx_x_R = f4.Histo2D(models[0], "k_vtx_x_R", "k_vtx_x_L")
h_vtx_y_L_vs_vtx_y_R = f4.Histo2D(models[1], "k_vtx_y_R", "k_vtx_y_L")

# Line 1231:
#            (k.th_x_L, k.vtx_x_L)
#            (k.th_x_R, k.vtx_x_R)
#            (k.th_y_L, k.vtx_y_L)
#            (k.th_y_R, k.vtx_y_R)
models = list(map(ROOT.ROOT.RDF.TH2DModel, [
    ("h_vtx_x_L_vs_th_x_L", ";#theta_{x}^{L};x^{*,L}", 100, -600E-6, +600E-6, 100, -0.5, +0.5),
    ("h_vtx_x_R_vs_th_x_R", ";#theta_{x}^{R};x^{*,R}", 100, -600E-6, +600E-6, 100, -0.5, +0.5),
    ("h_vtx_y_L_vs_th_y_L", ";#theta_{y}^{L};y^{*,L}", 100, -600E-6, +600E-6, 100, -0.5, +0.5),
    ("h_vtx_y_R_vs_th_y_R", ";#theta_{y}^{R};y^{*,R}", 100, -600E-6, +600E-6, 100, -0.5, +0.5)
]))
h_vtx_x_L_vs_th_x_L = f4.Histo2D(models[0], "k_th_x_L", "k_vtx_x_L")
h_vtx_x_R_vs_th_x_R = f4.Histo2D(models[1], "k_th_x_R", "k_vtx_x_R")
h_vtx_y_L_vs_th_y_L = f4.Histo2D(models[2], "k_th_y_L", "k_vtx_y_L")
h_vtx_y_R_vs_th_y_R = f4.Histo2D(models[3], "k_th_y_R", "k_vtx_y_R")

# Line 1236:
#           (k.vtx_x_R - k.vtx_x_L)
#           (k.vtx_y_R - k.vtx_y_L)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_vtx_x_diffLR", ";x^{*,R} - x^{*,L}", 100, -0.5, +0.5),
    ("h_vtx_y_diffLR", ";y^{*,R} - y^{*,L}", 100, -0.5, +0.5)
]))
h_vtx_x_diffLR = f4.Histo1D(models[0], "k_vtx_x_diffLR");
h_vtx_y_diffLR = f4.Histo1D(models[1], "k_vtx_y_diffLR");

# Line 1239:
#           (k.th_x, k.vtx_x_R - k.vtx_x_L)
#           (k.th_y, k.vtx_y_R - k.vtx_y_L)
models = list(map(ROOT.ROOT.RDF.TH1DModel, [
    ("h_vtx_x_diffLR", ";x^{*,R} - x^{*,L}", 100, -0.5, +0.5),
    ("h_vtx_y_diffLR", ";y^{*,R} - y^{*,L}", 100, -0.5, +0.5)
]))
h_vtx_x_diffLR_vs_th_x = f4.Histo1D(models[0], "k_th_x", "k_vtx_x_diffLR");
h_vtx_y_diffLR_vs_th_y = f4.Histo1D(models[1], "k_th_y", "k_vtx_y_diffLR");

# Line 1245:
#           (k.vtx_x_R, k.vtx_x_R - k.vtx_x_L)
#           (k.vtx_y_R, k.vtx_y_R - k.vtx_y_L)
models = list(map(ROOT.ROOT.RDF.TH2DModel,[
    ("h_vtx_x_diffLR_vs_vtx_x_R", ";x^{*,R};x^{*,R} - x^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5),
    ("h_vtx_y_diffLR_vs_vtx_y_R", ";y^{*,R};y^{*,R} - y^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5)
]))
h_vtx_x_diffLR_vs_vtx_x_R = f4.Histo2D(models[0], "k_vtx_x_R", "k_vtx_y_diffLR");
h_vtx_y_diffLR_vs_vtx_y_R = f4.Histo2D(models[1], "k_vtx_y_R", "k_vtx_y_diffLR");

# Define normalization and norm_corr colums
r6 = f4.Define("norm_corr",     "getNorm_corr( timestamp )" )         .Define("normalization", "getNormalization( norm_corr )")


# calculate acceptance divergence correction
r7 = r6.Define("correction", "CalculateAcceptanceCorrectionsRDF( th_y_sign, kinematics, anal )")         .Define("corr",       "correction.corr")         .Define("div_corr",   "correction.div_corr")         .Define("one",        "One()")

# Line 1406
for bi in binnings:
    bis = binning_setup[bi]

    model = ROOT.RDF.TH1DModel("h_t_Nev_before", ";|t|;events per bin", bis.N_bins, bis.bin_edges)
    bh_t_Nev_before[bi] = r7.Histo1D(model, "k_t", "one");

    model = ROOT.RDF.TH1DModel("h_t_before", ";|t|", bis.N_bins, bis.bin_edges)
    bh_t_before[bi] = r7.Histo1D(model, "k_t", "one");


# Line 1412
model = ROOT.RDF.TH2DModel("h_th_y_vs_th_x_before", ";#theta_{x};#theta_{y}", 150, -300E-6, +300E-6, 150, -150E-6, +150E-6)
h_th_y_vs_th_x_before = r7.Histo2D(model, "k_th_x", "k_th_y", "one");

# Line 1414
# Filter skip
f5 = r7.Filter("! correction.skip", "acceptance correction")

# Line 1429
for bi in binnings:
    bis = binning_setup[bi]

    model = ROOT.RDF.TH1DModel("h_t_Nev_after_no_corr", ";|t|;events per bin", bis.N_bins, bis.bin_edges)
    bh_t_Nev_after_no_corr[bi] = f5.Histo1D(model, "k_t", "one");

    model = ROOT.RDF.TH1DModel("h_t_after_no_corr", ";|t|", bis.N_bins, bis.bin_edges)
    bh_t_after_no_corr[bi] = f5.Histo1D(model, "k_t", "one");

    model = ROOT.RDF.TH1DModel("h_t_after", ";|t|", bis.N_bins, bis.bin_edges)
    bh_t_after[bi] = f5.Histo1D(model, "k_t", "corr");

# Line 1435
model = ROOT.RDF.TH2DModel("h_th_y_vs_th_x_after", ";#theta_{x};#theta_{y}", 150, -300E-6, +300E-6, 150, -150E-6, +150E-6);
h_th_y_vs_th_x_after = f5.Histo2D(model, "k_th_x", "k_th_y", "div_corr");

# Line 1435
model = ROOT.RDF.TH2DModel("h_th_vs_phi_after", ";#phi;#theta", 50, -M_PI, +M_PI, 50, 150E-6, 550E-6);
h_th_vs_phi_after = f5.Histo2D(model, "k_phi", "k_th", "div_corr");

# Line 1441
# apply normalization
model = ROOT.ROOT.RDF.TH1DModel("h_t_normalized", ";|t|",128, 0., 4.)
bh_t_normalized_ob_1_30_02 = f5.Define("corr_norm", "corr * normalization")                                 .Histo1D(model, "k_t", "corr_norm")

# Line 1445
model = ROOT.RDF.TH2DModel("h_th_y_vs_th_x_normalized", ";#theta_{x};#theta_{y}", 150, -600E-6, +600E-6, 150, -600E-6, +600E-6);
h_th_y_vs_th_x_normalized = f5.Define("div_corr_norm", "correction.div_corr * normalization")                                 .Histo2D(model, "k_th_x", "k_th_y", "div_corr_norm");


# normalize histograms
for bi in binnings:
    bh_t_before[bi].Scale(1., "width");
    bh_t_after_no_corr[bi].Scale(1., "width");
    bh_t_after[bi].Scale(1., "width");

# TODO: Write rest of files
outF = ROOT.TFile.Open("distributions_pyrdf_spark.root", "recreate");
h_timestamp_dgn.Write()
outF.Close()
