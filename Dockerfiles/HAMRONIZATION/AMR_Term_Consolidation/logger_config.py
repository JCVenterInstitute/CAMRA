import logging

# Configure logger
def setup_logger(name="term_consolidation", log_file="term_consolidation.log"):
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # Prevent duplicate logging
    if not logger.handlers:
        # Create formatter
        formatter = logging.Formatter("[{levelname} {asctime}] {message}", style="{")

        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger
