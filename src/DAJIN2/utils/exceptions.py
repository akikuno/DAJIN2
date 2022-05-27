class InputFileError(Exception):
    pretext = "Input File Error"

    def __init__(self, message, *args):
        if self.pretext:
            message = f"{self.pretext}: {message}"
        super().__init__(message, *args)

