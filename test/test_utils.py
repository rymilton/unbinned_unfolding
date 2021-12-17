import os
import ROOT
import json

def get_combination(parms, parms_name):
    combined_parms = []
    if len(parms_name) == 0:
        return [' ']
    cur_parm = parms_name[0]
    parms_name.pop(0)
    perv_comb = get_combination(parms, parms_name)
    for parm_val in parms[cur_parm]:
        for comb in perv_comb:
            new_comb = comb + " " + cur_parm + "=" + parm_val
            combined_parms.append(new_comb)
    return combined_parms

def get_unfold(f):
    unfold = f.Get("unfold")
    h_unfolded = unfold.Hunfold()
    #u = ROOT.TVector(h_unfolded.GetNbinsX())
    u = [0]*(h_unfolded.GetNbinsX())
    for i in range(h_unfolded.GetNbinsX()):
        u[i] = h_unfolded.GetBinContent(i+1)

    return u

def get_unfold_overflow(f):
    unfold = f.Get("unfold")
    h_unfolded = unfold.Hunfold()
    u = [0]*(h_unfolded.GetNbinsX()+2)
    for i in range(h_unfolded.GetNbinsX()+2):
        u[i] = h_unfolded.GetBinContent(i)

    return u

def get_unfold2D(f):
    unfold = f.Get("unfold")
    h_unfolded = unfold.Hunfold()
    #u = ROOT.TVector(h_unfolded.GetNbinsX()*h_unfolded.GetNbinsY())
    u = [0]*(h_unfolded.GetNbinsX()*h_unfolded.GetNbinsY())
    i = 0
    for x in range(h_unfolded.GetNbinsX()):
        for y in range(h_unfolded.GetNbinsY()):
            u[i] = h_unfolded.GetBinContent(x+1, y+1)
            i += 1
    return u

def get_unfold3D(f):
    unfold = f.Get("unfold")
    h_unfolded = unfold.Hunfold()
    #u = ROOT.TVector(h_unfolded.GetNbinsX()*h_unfolded.GetNbinsY()*h_unfolded.GetNbinsZ())
    u = [0]*(h_unfolded.GetNbinsX()*h_unfolded.GetNbinsY()*h_unfolded.GetNbinsZ())
    i = 0
    for x in range(h_unfolded.GetNbinsX()):
        for y in range(h_unfolded.GetNbinsY()):
            for z in range(h_unfolded.GetNbinsZ()):
                u[i] = h_unfolded.GetBinContent(x+1, y+1, z+1)
                i += 1
    return u

def get_uncertainty(f):
    unfold = f.Get("unfold")
    h_unfolded = unfold.Hunfold()
    u = ROOT.TVector(h_unfolded.GetNbinsX())
    u = [0]*(h_unfolded.GetNbinsX())
    for i in range(h_unfolded.GetNbinsX()):
        u[i] = h_unfolded.GetBinError(i+1)
    return u

def get_field(filename, field_to_compare = ['unfold']):
    global comparing_fields
    f =  ROOT.TFile.Open(filename,"READ")
    u = {}
    for field in field_to_compare:
        if field in comparing_fields:
            tmp_u = comparing_fields[field](f)
            u[field] = tmp_u
        else:
            print("[ERROR] Specified field not present in the file for comparison")
            exit(1)

    return u


def write_field(all_output, ref_file_name):
    with open(ref_file_name, "w") as outfile: 
        json.dump(all_output, outfile, indent=4)


def delete_files():
    os.system("rm -f RooUnfoldTest.root")
    os.system("rm -f RooUnfoldTest.ps")


def compare(all_output, ref_file_name, test_name, allowed_difference):
    data = {}
    with open(ref_file_name, "r") as f:
        data = json.load(f)

    if len(data) != len(all_output):
        print("Number of parms not same")
        return 1
    for parm, value in data.items():
        for field, num_list in value.items():
            for i in range(len(num_list)):
                if abs(num_list[i] - all_output[parm][field][i]) > allowed_difference:
                    print("Error in {} in field {} with parm {} at index {} as {} should be {}".format(test_name, field, parm, i, num_list[i], all_output[parm][field][i]))
                    return 1
    return 0

def perform_test(parms, ref_file_name, test_name, field_to_compare, allowed_difference = 1, is_combined = False):
    all_output = {}
    combined_parm = []
    if is_combined:
        combined_parm = parms
    else:
        parms_name = list(parms.keys())
        parms_name.sort()
        combined_parm = get_combination(parms, parms_name)

    delete_files()
    for single_parm in combined_parm:
        command_str = "../build/RooUnfoldTest " +  single_parm
        os.system(command_str)
        u = get_field("RooUnfoldTest.root", field_to_compare)
        all_output[single_parm] = u
        delete_files()

    if compare(all_output, ref_file_name, test_name, allowed_difference) == 1:
        print("Test failed")
        exit(1)


comparing_fields = {"unfold": get_unfold, 
                    "uncertainty":get_uncertainty, 
                    "unfold2D":get_unfold2D, 
                    "unfold3D":get_unfold3D,
                    "unfoldoverflow":get_unfold_overflow
                }
