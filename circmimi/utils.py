
def add_prefix(filename, prefix):
    filename = "{}.{}".format(prefix, filename)
    filename = filename.lstrip('.')
    return filename
