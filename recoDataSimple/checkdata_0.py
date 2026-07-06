import ROOT

file = input("Enter the run number: ")
filename = "recoDataSimple_" + file + "_xtalMerging.root"
df = ROOT.RDataFrame("simpleEvent", filename)

# check run number is always the same
values_to_check = ["Event.run", "Event.nuclear", "Event.nuclearRaw"]
for value in values_to_check:
    min_val = df.Min(value)
    max_val = df.Max(value)

    if min_val.GetValue() == max_val.GetValue():
        print("ok.")
    else:
        print(f"not ok for {value}.")

# check position of the goniometer does not change
values_to_check_1 = ["GonioPos.x", "GonioPos.y", "GonioPos.z"]
for value in values_to_check_1:
    min_val = df.Min(value)
    max_val = df.Max(value)

    if min_val.GetValue() == max_val.GetValue():
        print("ok.")
    else:
        print(f"not ok for {value}.")

# check single track and multi hits numbers
values_to_check_2 = ["SingleTrack", "MultiHit"]
for value in values_to_check_2:
    min_val = df.Min(value)
    max_val = df.Max(value)

    print(f"{value}: [min {min_val.GetValue()} : max {max_val.GetValue()}]")

mismatches = df.Filter("SingleTrack != MultiHit").Count()

print(mismatches.GetValue())









