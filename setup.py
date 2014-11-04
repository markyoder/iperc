import distutils.core as dcp

'''
module2 = dcp.Extension('spam', include_dirs=['/home/myoder/Documents/Source','/usr/local/include'], libraries=['/usr/libs'], sources = ['spammodule.cpp'])
dcp.setup (name = 'spam', version='1.0', description='spam extension module', ext_modules = [module2])
'''

#module1 = dcp.Extension('ipercPy', include_dirs=['/home/myoder/Documents/Source','/usr/local/include'], libraries=['/usr/libs'],  sources=['lattice.cpp', 'site.cpp', 'cluster.cpp'])
module1 = dcp.Extension('ipercPy', include_dirs=['/home/myoder/Documents/Source','/usr/local/include'], libraries=['/usr/libs'],  sources=['ipercPy.cpp', 'lattice.cpp', 'site.cpp', 'cluster.cpp', 'bond.cpp'])
dcp.setup (name = 'ipercPy', version='1.0', description='Invasion percolation model', ext_modules = [module1])
