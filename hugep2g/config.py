import os

project_dir = os.path.dirname(os.path.abspath(__file__))
tmp_work_dir = "/tmp"

genewise_path = "genewise"
wu_blast_path = os.path.join(project_dir, 'dep', 'wublast')
genblasta_path = os.path.join(project_dir, 'dep', 'genblasta', 'genblasta_v1.0.4_linux_x86_64')
genblasta_wu_script = os.path.join(project_dir, 'dep', 'genblasta', 'genblasta_wu.py')
