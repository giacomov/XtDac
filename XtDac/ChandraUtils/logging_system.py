import logging

logging.basicConfig()


def get_logger(name):
    logger = logging.getLogger(name)

    logger.setLevel(logging.INFO)

    return logger
