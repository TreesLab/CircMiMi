import os.path


def add_prefix(filename, prefix, auto_dot=False):
    if auto_dot:
        prefix_dir, prefix_base = os.path.split(prefix)

        filename = "{}.{}".format(prefix_base, filename)
        filename = filename.lstrip('.')
        filename = os.path.join(prefix_dir, filename)
    else:
        filename = prefix + filename

    return filename
