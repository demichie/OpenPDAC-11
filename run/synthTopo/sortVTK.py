import os

def find_and_create_links(root_folder):
    vtk_files = {}
    top_dir = os.getcwd()
    os.chdir(root_folder)    
    current_dir = os.getcwd()
    
    # Iterate over subfolders
    for subdir in sorted(os.listdir(current_dir)):
        subdir_path = os.path.join(current_dir, subdir)
        if os.path.isdir(subdir_path):
            # Iterate over files in subfolder
            for file in os.listdir(subdir_path):
                if file.endswith('.vtk'):
                    base_name = os.path.splitext(file)[0]
                    vtk_files.setdefault(base_name, []).append((subdir, os.path.join(subdir_path, file)))

    
    
    # Create symbolic links
    for base_name, files in vtk_files.items():
        index = 0
        for subdir, file_path in files:
            new_name = f"{base_name}_{index:03d}.vtk"
            new_link_path = os.path.join(current_dir, new_name)
            os.symlink(file_path, new_link_path)
            print(f"Created symlink: {new_link_path} -> {file_path}")
            index += 1

    os.chdir(top_dir)    


# Example usage
root_folder = './postProcessing/isoSurfaces/'
find_and_create_links(root_folder)

root_folder = './postProcessing/cuttingPlane/'
find_and_create_links(root_folder)

root_folder = './postProcessing/cutPlane/'
find_and_create_links(root_folder)

