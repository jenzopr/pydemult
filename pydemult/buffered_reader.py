CHAR_PLUS = ord(b'+')
CHAR_NEWLINE = ord(b'\n')

def buffered_blob(handle, bufsize):
    backlog = b''
    blob = b''

    while True:
        blob = handle.read(bufsize)
        if blob == b'':
            # Could be end of files
            yield backlog
            break

        if backlog != b'':
            blob = backlog + blob

        i = blob.rfind(b'\n@', 0)

        # Make sure no quality line was found
        if (blob[i-1] == CHAR_PLUS) and (blob[i-2] == CHAR_NEWLINE):
            i = blob.rfind(b'\n@', 0, i-2)

        backlog = blob[i+1:len(blob)]
        yield blob[0:i]
