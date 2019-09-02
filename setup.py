import setuptools


setuptools.setup(name='growai',
      version="1.0.0",
      url = "https://github.com/nostrumbiodiscovery/growai", 
      description='Coupling active learning approaches with generative models to automatically generate and rank 1M compounds.',
      author='Daniel Soler',
      author_email='daniel.soler@nostrumbiodiscovery.com',
      install_requires=["tqdm", ],
      packages=setuptools.find_packages(),
      classifiers=[
       "Programming Language :: Python :: 3",
       "License :: OSI Approved :: MIT License",
       "Operating System :: OS Independent",
       ],
     )
