from setuptools import setup

with open("README.md", "r") as fh:
      long_description = fh.read()

setup(name='sciviewer',
      version='0.12',
      description='Tool for interactive exploration of 2D embeddings in Jupyter',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/colabobio/sciviewer',
      author='Andres Colubri, Dylan Kotliar',
      author_email='andres@broadinstitute.org, dylan_kotliar@hms.harvard.edu',
      license='MIT',
      packages=['sciviewer'],
      install_requires=[
          'py5',
      ],
      zip_safe=False)
