import re
import os

from mzlib.index import MemoryIndex

from .base import _PlainTextSpectralLibraryBackendBase


class TextSpectralLibrary(_PlainTextSpectralLibraryBackendBase):

    def _build_index(self):
        """
        Populate the spectrum index

        Returns
        -------
        n_spectra: int
            The number of entries read
        """

        #### Check that the spectrum library filename isvalid
        filename = self.filename

        #### Determine the filesize
        file_size = os.path.getsize(filename)

        with open(filename, 'r') as infile:
            state = 'header'
            spectrum_buffer = []
            n_spectra = 0
            start_index = 0
            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''

            # Required for counting file_offset manually (LF vs CRLF)
            infile.readline()
            file_offset_line_ending = len(infile.newlines) - 1
            infile.seek(0)
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                # file_offset = infile.tell()
                file_offset += len(line) + file_offset_line_ending

                line = line.rstrip()
                if state == 'header':
                    if re.match('MS:1003061|spectrum name=', line):
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                    else:
                        continue
                if state == 'body':
                    if len(line) == 0:
                        continue
                    if re.match('MS:1003061|spectrum name=', line):
                        if len(spectrum_buffer) > 0:
                            self.index.add(
                                number=n_spectra + start_index,
                                offset=spectrum_file_offset,
                                name=spectrum_name,
                                analyte=None)
                            n_spectra += 1
                            spectrum_buffer = []
                            #### Commit every now and then
                            if n_spectra % 1000 == 0:
                                self.index.commit()
                            #     percent_done = int(
                            #         file_offset/file_size*100+0.5)
                            #     eprint(str(percent_done)+"%..",
                            #            end='', flush=True)

                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match('MS:1003061|spectrum name=(.+)', line).group(1)

                    spectrum_buffer.append(line)

            self.index.add(
                number=n_spectra + start_index,
                offset=spectrum_file_offset,
                name=spectrum_name,
                analyte=None)
            self.index.commit()
            n_spectra += 1

            #### Flush the index
            self.index.commit()

        return n_spectra

    def _get_lines_for(self, offset):
        with open(self.filename, 'r') as infile:
            infile.seek(offset)
            state = 'body'
            spectrum_buffer = []

            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''
            for line in infile:
                line_beginning_file_offset = file_offset
                file_offset += len(line)
                line = line.rstrip()
                if state == 'body':
                    if len(line) == 0:
                        continue
                    if re.match('MS:1003061|spectrum name=', line):
                        if len(spectrum_buffer) > 0:
                            return spectrum_buffer
                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match(
                            'MS:1003061|spectrum name=(.+)', line).group(1)
                    spectrum_buffer.append(line)

            #### We will end up here if this is the last spectrum in the file
            return spectrum_buffer
