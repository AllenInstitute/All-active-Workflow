from setuptools import setup,find_packages

setup(name='ateamopt',
      version='0.1',
      description='Ateam all-active optimization toolbox',
      author='Ani Nandi',
      author_email='anin@alleninstitute.org',
      packages=find_packages(),
      scripts=['ateamopt/jobscript/submit_opt_jobs'],
      entry_points={
        'console_scripts':[
            'launch_optimjob = ateamopt.jobscript.launch_optimjob:main'
        ]
      },
      platforms='any'
    )

