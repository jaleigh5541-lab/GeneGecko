from io import BytesIO


class FileWrapper:
    """Adapts a filename + raw bytes to the .name/.read()/.seek() interface
    expected by Prot_modules and intron_modules."""

    def __init__(self, filename: str, data: bytes) -> None:
        self.name = filename
        self._stream = BytesIO(data)

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)

    def seek(self, pos: int, whence: int = 0) -> int:
        return self._stream.seek(pos, whence)

    def tell(self) -> int:
        return self._stream.tell()

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True
