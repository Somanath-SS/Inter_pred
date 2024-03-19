import pathlib as pl
import os

class PathHandler():
    def __init__(self, folder_path):
        self.path = pl.Path(folder_path)

    def get_pdb_files(self):
        pdb_files = list(self.path.glob('**/*.pdb'))
        return pdb_files, len(pdb_files)
    
if __name__ == "__main__":
    handler = PathHandler(r'C:\Users\Damian\PycharmProjects\Inter_pred\sample_data')
    pdb_files = handler.get_pdb_files()
    for pdb_file in pdb_files:
        print(pdb_file.absolute)