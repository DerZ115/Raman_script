import os
import re

data_path = "C:/Users/Daniel/Desktop/Raman_script/spectra/Behandlungsversuche_neu/E2_txt"
group_name = "E2"

files = os.listdir(data_path)
n = len(files)
digits = len(str(n))

files.sort(key=lambda x: int(re.findall(r'\d+', x)[-1]))

files_new = []

for i in range(n):
    new_filename = group_name + "_" + str(i+1).zfill(digits) + ".TXT"
    # os.rename(os.path.join(data_path, files[i]), os.path.join(data_path, new_filename))
    print(f"{i}: {files[i]} -> {new_filename}")
    




