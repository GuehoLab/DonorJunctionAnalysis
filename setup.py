# setup.py

from setuptools import setup, find_packages

setup(
    name='donor-junction-analysis',
    version='0.1.0',
    description='A bioinformatics pipeline for analyzing mutations and donor junctions.',
    author='Your Name',
    
    # List of requirements (optional, as they are mostly in environment.yml)
    install_requires=[
        'pysam==0.23.3',
        'matplotlib==3.10.7',
        'edlib==1.3.9',
    ],


    # Automatically find the 'src' directory as the package root
    package_dir={'': 'scr'},
    packages=find_packages(where='scr'),
    
    # entry_points에서는 이제 'src' prefix 없이 모듈명을 사용합니다.
    entry_points={
        'console_scripts': [
            'donor_junction_analysis = DonorJunctionAnalysis:main',
            'junction_trim_insert = junction_trim_insert:main',
            'junction_analyze_mut = junction_analyze_mut:main'
        ],
    },
)