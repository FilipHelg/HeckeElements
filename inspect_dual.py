import struct
path='dual_kl_bulk_n5.bin'
try:
    with open(path,'rb') as f:
        header=f.read(32)
        if len(header)<32:
            print('Header too short')
            raise SystemExit
        magic=header[0:8].decode('ascii',errors='replace')
        version, n, count, involutionCount, recordFormat = struct.unpack('<5I', header[8:28])
        print('magic=',magic,'version=',version,'n=',n,'count=',count,'involutionCount=',involutionCount,'recordFormat=',recordFormat)
        found=False
        for rec in range(involutionCount):
            raw=f.read(8)
            if len(raw)<8:
                print('truncated at rec',rec)
                break
            x, supportCount = struct.unpack('<II', raw)
            if x==7:
                found=True
                print('Found record for x=7 with supportCount=',supportCount)
                # read support entries and print first few
                for s in range(supportCount):
                    raw=f.read(5)
                    if len(raw)<5:
                        print('truncated support entry')
                        break
                    z, termCount = struct.unpack('<IB', raw)
                    print(' support z=',z,'termCount=',termCount)
                    terms=[]
                    for t in range(termCount):
                        rawt=f.read(5)
                        if len(rawt)<5:
                            break
                        exponent, coeff = struct.unpack('<bi', rawt)
                        terms.append((exponent,coeff))
                    print('  terms=',terms[:10])
                break
            else:
                # skip supportCount entries
                for s in range(supportCount):
                    raw=f.read(5)
                    if len(raw)<5:
                        break
                    z, termCount = struct.unpack('<IB', raw)
                    f.read(5*termCount)
        if not found:
            print('No record for x=7 in dual cache')
except FileNotFoundError:
    print('dual cache not found')
except Exception as e:
    print('Error:',e)
