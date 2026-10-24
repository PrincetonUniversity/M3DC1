
import numpy as np
import sys
import matplotlib.pyplot as plt

def chunkstring(string, length):
  return (string[0+i:length+i] for i in range(0, len(string), length))

def chunklist(string, length):
  return [item for item in chunkstring(string, length) if item!='\n']

def readg(fname):
  f = open(fname, 'r')

  g={}
  data_order=['fpol','pres','ffprim','pprim','psirz','qpsi','nblm','rzbbbs','rzlim']
  j=0
  k=0
  l=0

  for i, line in enumerate(f):

    if i==0:
      header = line.split()
      g['header']=header
      g['nw']=int(header[-2])
      g['nh']=int(header[-1])
      nw=g['nw']
      nh=g['nh']
      #print(header)
      for item in data_order:
        if item=='psirz':
          g[item]=np.zeros((nh,nw))
        elif item in ['nbbbs','limitr','rzbbbs','rzlim']: # we dont know this size yet
          pass
        else:
          g[item]=np.zeros(nw) # all the 1D profs
    elif i==1:
      list_of_chunks=chunklist(line, 16)
      #print(list_of_chunks)
      #for item in list_of_chunks:
      #  print(item)
      g['rdim']=float(list_of_chunks[0])
      g['zdim']=float(list_of_chunks[1])
      g['rcentr']=float(list_of_chunks[2])
      g['rleft']=float(list_of_chunks[3])
      g['zmid']=float(list_of_chunks[4])

      rleft=g['rleft']
      rdim=g['rdim']
      zdim=g['zdim']

      g['rgrid']=np.linspace(rleft, rleft+rdim, nw)
      g['zgrid']=np.linspace(-zdim/2, zdim/2, nh)

    elif i==2:
      list_of_chunks=chunklist(line, 16)
      g['rmaxis']=float(list_of_chunks[0])
      g['zmaxis']=float(list_of_chunks[1])
      g['simag']=float(list_of_chunks[2])
      g['sibry']=float(list_of_chunks[3])
      g['bcentr']=float(list_of_chunks[4])
    elif i==3:
      list_of_chunks=chunklist(line, 16)
      g['current']=float(list_of_chunks[0])
      g['simag']=float(list_of_chunks[1])
      g['xdum']=float(list_of_chunks[2])
      g['rmaxis']=float(list_of_chunks[3])
      g['xdum']=float(list_of_chunks[4])
    elif i==4:
      list_of_chunks=chunklist(line, 16)
      g['zmaxis']=float(list_of_chunks[0])
      g['xdum']=float(list_of_chunks[1])
      g['sibry']=float(list_of_chunks[2])
      g['xdum']=float(list_of_chunks[3])
      g['xdum']=float(list_of_chunks[4])

 # rmaxis=g(1);zmaxis=g(2);simag=g(3);sibry=g(4);bcentr=g(5);

    else:
      for chunk in chunklist(line, 16):
        if l<len(data_order):
          #print(data_order[l],j,k,l,chunk)
          if data_order[l]=='psirz':
            g[data_order[l]][k][j]=float(chunk)
            j=j+1
            if j==nw:
              j=0
              k=k+1
            if k==nh:
              j=0
              k=0
              l=l+1
          elif data_order[l]=='nblm':
            line_items=line.split()
            g['nbbbs']=int(line_items[0])
            g['rbbbs']=np.zeros(g['nbbbs'])
            g['zbbbs']=np.zeros(g['nbbbs'])
            g['limitr']=int(line_items[1])
            g['rlim']=np.zeros(g['limitr'])
            g['zlim']=np.zeros(g['limitr'])
            j=nw # this moves on to next var
          elif data_order[l]=='rzbbbs':
            if k==0:
              g['rbbbs'][j]=float(chunk)
              k=1
            else:
              g['zbbbs'][j]=float(chunk)
              k=0
              j=j+1
            if j==g['nbbbs']:
              j=0
              k=0
              l=l+1
          elif data_order[l]=='rzlim':
            if k==0:
              g['rlim'][j]=float(chunk)
              k=1
            else:
              g['zlim'][j]=float(chunk)
              k=0
              j=j+1
            if j==g['limitr']:
              j=0
              k=0
              l=l+1
          else: # all the 1D profiles
            g[data_order[l]][j]=float(chunk)
            j=j+1
          if j==nw and data_order[l] not in ['nbbbs','limitr','rzbbbs','rzlim','psirz']:
            j=0
            k=0
            l=l+1
  f.close()

  return g

if __name__ == '__main__':
  print("file:", sys.argv[1])
  g = readg(sys.argv[1])
  # then we do something with g, like plot parts like readg.m
