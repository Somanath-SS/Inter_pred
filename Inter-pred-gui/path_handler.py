import pathlib as pl
import os

class PathHandler():
    """
    This class is used to handle the path functionality
    
    Args:
        folder_path (str): The path to the folder
        
    Functions:
        get_pdb_files: This function is used to get the pdb files
    
    """
    def __init__(self, folder_path):
        self.path = pl.Path(folder_path)

    def get_pdb_files(self):
        """
        This function is used to get the pdb files
        
        Args:
            None: This function does not take any arguments
        
        Returns:
            pdb_files (list): The list of pdb files
            len(pdb_files) (int): The number of pdb files
        
        Raises:
            None: This function does not raise any exceptions
        """
        pdb_files = list(self.path.glob('**/*.pdb'))
        return pdb_files, len(pdb_files)
    
if __name__ == "__main__":
    handler = PathHandler(r'C:\Users\Damian\PycharmProjects\Inter_pred\sample_data')
    pdb_files = handler.get_pdb_files()
    for pdb_file in pdb_files:
        print(pdb_file.absolute)