#%%
import os
import zipfile
import shutil
import subprocess
import requests
from sys import platform

group = input("Group number ? >> ");
login1 = input("Login of the first student ? >> ")
login2 = input("Login of the second student ? >> ")
lsgroup = "%03d" % int(group)
login1 = login1.lower()
login2 = login2.lower()
idname = "group" + lsgroup + '-' + login1 + '-' + login2
print("Your file name must be : " + idname)

# Configuration A MODIFIER PAR L'ETUDIANT !!!
# =============================================
# project_path = "group001-vlegat-jfremacle.zip"
# rapport_path = "group001-vlegat-jfremacle-rapport.zip"
project_path = idname + ".zip"
rapport_path = idname + "-rapport.zip"
Windows_TDM = (platform == "win32")  # True if using Windows + TDM-GCC, False otherwise
# =============================================


all_illegal = [
    "/gmsh/",
    "/build/",
    "/.",
]

project_illegal = [
    "/glfw/",
    "CMakeCache.txt",
    "Makefile",
    "CMakeFiles",
    "ProjectPostProcessor/",
    "ProjectPreProcessor/",
    "/data/"
]

project_allow = [
    ".c",
    ".cpp",
    ".h",
    ".hpp",
    ".txt",
    ".md"
]

# %% File name validation
project_fname = os.path.basename(project_path)
rapport_fname = os.path.basename(rapport_path)
assert project_fname.endswith(".zip"), "Project file must be a zip file"
assert rapport_fname.endswith(".zip"), "Rapport file must be a zip file"
project_parts = project_fname[:-4].split("-")
rapport_parts = rapport_fname[:-4].split("-")
assert len(project_parts) == 3, "Project file name must be groupXXX-YYY-ZZZ.zip"
assert len(rapport_parts) == 4, "Rapport file name must be groupXXX-YYY-ZZZ-rapport.zip"
assert project_parts[0] == rapport_parts[0], "file names must start with 'groupXXX'"
assert project_parts[0][5:] == rapport_parts[0][5:] == "%0.3d" % int(project_parts[0][5:]), "group number must be a 3 digit number"
assert project_parts[1] == rapport_parts[1], "group names must be the same"
assert project_parts[2] == rapport_parts[2], "group names must be the same"
assert rapport_parts[3] == "rapport", "rapport file name must end with '-rapport'"
assert project_parts[1] != "vlegat", "Non, c'est votre identifiant ucl qu'il faut mettre, vlegat-jfremacle c'est pour l'exemple"
assert project_parts[2] != "jfremacle", "Non, c'est votre identifiant ucl qu'il faut mettre, vlegat-jfremacle c'est pour l'exemple"

# %% Rapport content validation
if os.path.exists('./tmp-rapport'):
    shutil.rmtree('./tmp-rapport')

illegal = 0
penalty = 0
os.makedirs('./tmp-rapport')
with zipfile.ZipFile(rapport_path, 'r') as zip:
    for file in zip.namelist():
        nstrip = 0
        if file.startswith(project_fname[:-4]):
            nstrip = len(project_fname[:-4]) + 1
        if any(i in file for i in all_illegal):
            print(f"WARNING ERROR : illegal file : {file}")
            illegal += 1; penalty += 10
            continue
        print(file)
        if file.endswith('/'):
            os.makedirs('./tmp-rapport/' + file[nstrip:], exist_ok=True)
            continue
        with zip.open(file) as zf, open('./tmp-rapport/' + file[nstrip:], 'wb') as f:
            shutil.copyfileobj(zf, f)

if not os.path.exists('./Validated'):
    os.makedirs('./Validated')

with zipfile.ZipFile(f"./Validated/{rapport_fname}", 'w', zipfile.ZIP_DEFLATED) as new_zip:
    print("Creating new zip file with files")
    for root, dirs, files in os.walk('./tmp-rapport'):
        for file in files:
            new_zip.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), "./tmp-rapport"))

# %% Project content validation
if os.path.exists('./tmp'):
    shutil.rmtree('./tmp')

os.makedirs('./tmp')
with zipfile.ZipFile(project_path, 'r') as zip:
    for file in zip.namelist():
        nstrip = 0
        if file.startswith(project_fname[:-4]):
            nstrip = len(project_fname[:-4]) + 1
        if not (file.startswith('Project/') or file.startswith(project_fname[:-4])):
            # print(f"WARNING ERROR : illegal file : {file}")
            illegal += 1; penalty += 1
            continue
        if not any(file.endswith(ext) for ext in project_allow) and '.' in file:
            # print(f"WARNING ERROR : illegal file : {file}")
            illegal += 1; penalty += 1
            continue
        if any(i in file for i in all_illegal):
            # print(f"WARNING ERROR : illegal file : {file}")
            illegal += 1; penalty += 10
            continue
        if any(i in file for i in project_illegal):
            # print(f"WARNING ERROR : illegal file : {file}")
            illegal += 1; penalty += 1
            continue
        if file.endswith('/'):
            os.makedirs('./tmp/' + file[nstrip:], exist_ok=True)
            continue
        with zip.open(file) as zf, open('./tmp/' + file[nstrip:], 'wb') as f:
            shutil.copyfileobj(zf, f)

# %% Project compilation
if os.path.exists('./tmp-compil'):
    shutil.rmtree('./tmp-compil')
shutil.copytree('./tmp', './tmp-compil')

if Windows_TDM:
    exec = "myFem.exe"
    res = subprocess.run("cd ./tmp-compil/Project && mkdir build && cd build && cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -G \"MinGW Makefiles\" && cmake --build .", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
else:
    exec = "./myFem"
    res = subprocess.run("cd ./tmp-compil/Project && mkdir build && cd build && cmake .. && cmake --build .", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

print(res.stdout.decode('utf-8'))
print(res.stderr.decode('utf-8'))
if not os.path.exists(f"./tmp-compil/Project/build/{exec}"):
    raise Exception("Compilation failed")
print("Compilation successful")

# %% Project execution
os.makedirs("./tmp-compil/Project/data", exist_ok=True)
r1 = requests.get("https://raw.githubusercontent.com/MiguelDLC/femtools/main/data/mesh.txt")
r2 = requests.get("https://raw.githubusercontent.com/MiguelDLC/femtools/main/data/problem.txt")
if r1.status_code != 200:
    raise Exception("Failed to download mesh.txt")
if r2.status_code != 200:
    raise Exception("Failed to download problem.txt")
with open("./tmp-compil/Project/data/mesh.txt", "wb") as f:
    f.write(r1.content)
with open("./tmp-compil/Project/data/problem.txt", "wb") as f:
    f.write(r2.content)

res = subprocess.run(f"cd ./tmp-compil/Project/build && {exec}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
print(res.stdout.decode('utf-8'))
print(res.stderr.decode('utf-8'))
if not os.path.exists("./tmp-compil/Project/data/UV.txt"):
    raise Exception("Execution failed")
print("Execution successful")

# %%
if not os.path.exists('./Validated'):
    os.makedirs('./Validated')
with zipfile.ZipFile(f"./Validated/{project_fname}", 'w', zipfile.ZIP_DEFLATED) as new_zip:
    print("Creating new zip file with files")
    for root, dirs, files in os.walk('./tmp'):
        for file in files:
            print("   ",  os.path.relpath(os.path.join(root, file), "./tmp"))
            new_zip.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), "./tmp"))

print("\n========== Validation summary ==========")
print(f"    Compilation and execution successful")
if illegal == 0:
    print(f"    No illegal files found, you are good to go!")
else:
    print(f"    Total illegal files : {illegal} : {-penalty/10}/30 points")
    print("    Please check the files above and remove them before submitting again.")
    print("    Minor things like a readme are ok. Light but unnecessary files can be forgiven, but please remove them if possible.")
    print("    Major things like a build folder or gmsh folder will be heavily penalized.")
print("Clean files written to ./Validated") 

# %%