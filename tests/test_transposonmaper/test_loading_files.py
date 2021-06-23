import os
from satay.transposonmapping.importing.loading_files import access_files


def test_access_files():
    """Test if the function acces_file is giving a correct file according the extension provided (positive control)"""

    output = access_files(file_path=None, extension=".gff3")
    path = os.path.join("satay/data_files/")
    filename = "Saccharomyces_cerevisiae.R64-1-1.99.gff3"
    file = os.path.join(path, filename)

    assert output == file, "Should be %s" % file

