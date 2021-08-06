from collections import namedtuple

Species = namedtuple("Species", ['key', 'name', 'fullname', 'tax_id'])

species_list = {
    'ath': Species('ath', 'arabidopsis_thaliana', 'Arabidopsis thaliana', '3702'),
    'bmo': Species('bmo', 'bombyx_mori', 'Bombyx mori', '7091'),
    'bta': Species('bta', 'bos_taurus', 'Bos taurus', '9913'),
    'cel': Species('cel', 'caenorhabditis_elegans', 'Caenorhabditis elegans', '6239'),
    'cfa': Species('cfa', 'canis_familiaris', 'Canis familiaris', '9615'),
    'cgr': Species('cgr', 'cricetulus_griseus_chok1gshd', 'Cricetulus griseus', '10029'),
    'dre': Species('dre', 'danio_rerio', 'Danio rerio', '7955'),
    'dme': Species('dme', 'drosophila_melanogaster', 'Drosophila melanogaster', '7227'),
    'gga': Species('gga', 'gallus_gallus', 'Gallus gallus', '9031'),
    'hsa': Species('hsa', 'homo_sapiens', 'Homo sapiens', '9606'),
    'mmu': Species('mmu', 'mus_musculus', 'Mus musculus', '10090'),
    'osa': Species('osa', 'oryza_sativa', 'Oryza sativa', '4530'),
    'ola': Species('ola', 'oryzias_latipes', 'Oryzias latipes', '8090'),
    'oar': Species('oar', 'ovis_aries', 'Ovis aries', '9940'),
    'rno': Species('rno', 'rattus_norvegicus', 'Rattus norvegicus', '10116'),
    'ssc': Species('ssc', 'sus_scrofa', 'Sus scrofa', '9823'),
    'tgu': Species('tgu', 'taeniopygia_guttata', 'Taeniopygia guttata', '59729'),
    'xtr': Species('xtr', 'xenopus_tropicalis', 'Xenopus tropicalis', '8364')
}
