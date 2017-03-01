import os

cpp_extensions = (".cpp", ".cxx", ".c++", ".cc", ".cp", ".c", ".C", ".h", ".h++", ".hpp", ".hxx", ".hh")

for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(cpp_extensions):
            print(file)
            os.system("clang-format -i -style=file " + root + "/" + file)