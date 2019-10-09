from collections import namedtuple

Species = namedtuple("Species", ['key', 'name', 'fullname'])

species_list = {
    'ath': Species('ath', 'arabidopsis_thaliana', 'Arabidopsis thaliana'),
    'bmo': Species('bmo', 'bombyx_mori', 'Bombyx mori'),
    'bta': Species('bta', 'bos_taurus', 'Bos taurus'),
    'cel': Species('cel', 'caenorhabditis_elegans', 'Caenorhabditis elegans'),
    'cfa': Species('cfa', 'canis_familiaris', 'Canis familiaris'),
    'dre': Species('dre', 'danio_rerio', 'Danio rerio'),
    'dme': Species('dme', 'drosophila_melanogaster', 'Drosophila melanogaster'),
    'gga': Species('gga', 'gallus_gallus', 'Gallus gallus'),
    'hsa': Species('hsa', 'homo_sapiens', 'Homo sapiens'),
    'mmu': Species('mmu', 'mus_musculus', 'Mus musculus'),
    'osa': Species('osa', 'oryza_sativa', 'Oryza sativa'),
    'ola': Species('ola', 'oryzias_latipes', 'Oryzias latipes'),
    'oar': Species('oar', 'ovis_aries', 'Ovis aries'),
    'rno': Species('rno', 'rattus_norvegicus', 'Rattus norvegicus'),
    'ssc': Species('ssc', 'sus_scrofa', 'Sus scrofa'),
    'tgu': Species('tgu', 'taeniopygia_guttata', 'Taeniopygia guttata'),
    'xtr': Species('xtr', 'xenopus_tropicalis', 'Xenopus tropicalis')
}
