from setuptools import setup,find_packages

setup(
        name='covidquant',    # This is the name of your PyPI-package.
        version='1.0.0',    # Update the version number for new releases
        scripts=['CovidQuant.py',], # The name of your script, and also the command you'll be using for calling it
        description='CovidQuant: Quantify the abundance of different lineages in a mixed SARS-CoV-2 sequencing dataset.',
        long_description='CovidQuant: Quantify the abundance of different lineages in a mixed SARS-CoV-2 sequencing dataset.',
        long_description_content_type='text/plain',
        url='https://github.com/niemasd/CovidQuant',
        author='Niema Moshiri',
        author_email='niemamoshiri@gmail.com',
        packages=find_packages(),
        zip_safe = False,
        install_requires=['sklearn', 'pysam'],
        include_package_data=True
)
