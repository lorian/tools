import os



def file_split( filename, split_size, split_string ):
    """ Filename to split,
        splits when size is greater than split_size
        splits only on lines starting with split_string """

    f = open( filename, 'rt')

    file_count = 0
    line_count = 0
    output_fname_string = "{0}_{1}"

    output = open( output_fname_string.format( filename, file_count ), 'wt')
    for line in f:
        if line.startswith( split_string ):
            output.flush()
            fs = os.fstat( output.fileno() )
            if fs.st_size > split_size:
                output.close()
                file_count += 1
                output = open( output_fname_string.format( filename, file_count ), 'wt')

        output.write( line )
        line_count += 1
        if line_count % 10000 == 0:
            print "Lines processed: {0}".format( line_count )

    output.close()




