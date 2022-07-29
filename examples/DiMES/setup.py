import os


os.system("python smalllargedotsDIMES.py")
os.system("python setup_geom_DIMES.py")
os.system("python particles_distribution_generator.py")
os.system("python input.py")
os.system("./GITR_no_geom_hash")