import os

data_path = "C:/Users/Daniel/Desktop/Raman_script/spectra/HL_Behandlungsversuche/Resveratrol2/"
group_name = "R2"

files = os.listdir(data_path)
n = len(files)
digits = len(str(n))

files.sort(key=lambda x: int(x[8:-4]))

files_new = []

for i in range(n):
    new_filename = group_name + "_" + str(i+1).zfill(digits) + ".TXT"
    os.rename(data_path + files[i], data_path + new_filename)
    print(f"{i}: {new_filename}")
    




