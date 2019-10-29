#Borja :)
import sys, os, subprocess

def executecmd(args, out = None):
    '''
    :param path: LIST of args
    :return:
    '''
    if out is None:
        subprocess.call(args)
    else:
        with open(out, 'w+') as fpout:
            subprocess.call(args, stdout = fpout)

def __buildIndex(file):
    folder = '/'.join(file.split('/')[:-1])
    cmd = ['sga','index','-p',folder,file]
    print('Command: ',cmd)
    executecmd(cmd)

def __rmDup(file):
    extra = len(file.split('.')[-1])
    folder = '/'.join(file.split('/')[:-1])
    cmd = ['sga','rmdup','-p',folder,'-o',file[:-extra],file]
    print('Command: ',cmd)
    executecmd(cmd)

if __name__=='__main__':
    print('Unitig duplicates removal')

    assert len(sys.argv) > 1

    unitigsFile = sys.argv[1]
    print('Processing: ',unitigsFile)
    __buildIndex(unitigsFile)
    __rmDup(unitigsFile)
