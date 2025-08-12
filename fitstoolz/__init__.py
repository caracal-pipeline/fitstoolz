import logging

def get_logger(name, level="INFO"):

    if isinstance(level, str):
        level = getattr(logging, level, 10)

    format_string = '%(asctime)s-%(name)s-%(levelname)-8s| %(message)s'
    # set up logging to file - see previous section for more details
    logging.basicConfig(level=level,
                        format=format_string,
                        datefmt='%m:%d %H:%M:%S')

    return logging.getLogger(name)


LOG = get_logger("fitstoolz")