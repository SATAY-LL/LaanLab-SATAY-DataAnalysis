import os
import pkg_resources

from satay.transposonmapping.properties.get_gene_position import gene_position


def test_packaged_data():
    """Test default arguments """

    data_path = pkg_resources.resource_filename("satay", "data_files/")
    data_file = "Saccharomyces_cerevisiae.R64-1-1.99.gff3"
    filepath = os.path.join(data_path, data_file)

    gene_dict = gene_position(gff_file=filepath, get_dict=True)

    assert isinstance(gene_dict, dict), "Expected dict type"
    assert len(gene_dict) == 6600, "Expected dictionary of size 6600"
    assert "YBR247C" in gene_dict, "Expected YBR247C key to be present"
