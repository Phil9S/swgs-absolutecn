import ruamel.yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString, DoubleQuotedScalarString

yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml.explicit_start = True

def get_input(x,y):
    new_input = str(input(f'Set {y} (current: {x}): '))
    if not new_input:
        return x
    else:
        return new_input

# Config file updates
print('Update config.yaml values (Return no value to keep the current value)')
file_name = 'config/config.yaml'
config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))

# Update config values
config['samplesheet'] = DoubleQuotedScalarString(get_input(config['samplesheet'],'samplesheet'))
config['bins'] = int(get_input(config['bins'],'bins'))
config['out_dir'] = DoubleQuotedScalarString(get_input(config['out_dir'],'out_dir'))
config['project_name'] = DoubleQuotedScalarString(get_input(config['project_name'],'project_name'))

# Save updated yaml
with open('config/config.yaml', 'w') as fp:
    yaml.dump(config, fp)

# Slurm profile updates
print('Update slurm profile config.yaml values (Return no value to keep the current value)')
file_name = 'profile/slurm/config.yaml'
config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))

# Update jobs
config['jobs'] = int(get_input(config['jobs'],'jobs'))

# Update slurm partition and acc
clustervals = config['cluster'].split(' ')
clustervals[2] = get_input(clustervals[2],'account')
clustervals[4] = get_input(clustervals[4],'partition')
config['cluster'] = DoubleQuotedScalarString(' '.join([str(elem) for elem in clustervals]))

# Update slurm default resources
cpu_res = config['default-resources'][0].split('=')
cpu_res[1] = get_input(cpu_res[1],'default cpus per job')
config['default-resources'][0] = '='.join([str(elem) for elem in cpu_res])

mem_res = config['default-resources'][1].split('=')
mem_res[1] = get_input(mem_res[1],'default mem (mb) per job')
config['default-resources'][1] = '='.join([str(elem) for elem in mem_res])

time_res = config['default-resources'][2].split('=')
time_res[1] = get_input(time_res[1],'default time (mins) per job')
config['default-resources'][2] = '='.join([str(elem) for elem in time_res])

# Save updated yaml
with open('profile/slurm/config.yaml', 'w') as fp:
    yaml.dump(config, fp)
