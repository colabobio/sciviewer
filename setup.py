from setuptools import setup

with open("README.md", "r") as fh:
      long_description = fh.read()

setup(name='umap_explorer',
      version='0.1',
      description='Tool for interactive exploration of 2D embeddings in Jupyter',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/broadinstitute/embedview',
      author='Andres Colubri, Dylan Kotliar',
      author_email='andres@broadinstitute.org, dylan_kotliar@hms.harvard.edu',
      license='MIT',
      packages=['umap_explorer'],
      install_requires=[
          'py5',
      ],
      zip_safe=False)
