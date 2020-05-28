import os.path


def add_prefix(filename, prefix):
    prefix_dir, prefix_base = os.path.split(prefix)

    filename = "{}.{}".format(prefix_base, filename)
    filename = filename.lstrip('.')
    filename = os.path.join(prefix_dir, filename)

    return filename
