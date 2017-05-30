import argparse

class TempyParser(object):

    
    def __init__(self,
                 args=None):
        self.args = args

        
    def generate_args(self):
        parser = argparse.ArgumentParser()
        #
        job_title = parser.add_argument_group()
        job_title.add_argument(
            '-j',
            '--job_title',
            help='Short description of job',
            metavar='Job title',
            type=str,
            dest='job_title',
            default=None,
            required=False)
       #
        inp_map = parser.add_argument_group()
        inp_map.add_argument(
            '-m',
            '--map',
            help=('Input map'),
            metavar='Input map',
            type=str,
            dest='inp_map',
            default=None,
            required=False)
        inp_map.add_argument(
            '-m1',
            '--map1',
            help=('Input map1'),
            metavar='Input map1',
            type=str,
            dest='inp_map1',
            default=None,
            required=False)
        inp_map.add_argument(
            '-m2',
            '--map2',
            help=('Input map2'),
            metavar='Input map2',
            type=str,
            dest='inp_map2',
            default=None,
            required=False)

        search_model = parser.add_argument_group()
        search_model.add_argument(
	    '-p',
            '--pdb',
            help=('Input model'),
            type=str,
            default=None,
	    dest='pdb',
	    required=False)
        search_model.add_argument(
        '-p1',
            '--pdb1',
            help=('Input model1'),
            type=str,
            default=None,
        dest='pdb1',
        required=False)
        search_model.add_argument(
        '-p2',
            '--pdb2',
            help=('Input model2'),
            type=str,
            default=None,
        dest='pdb2',
        required=False)
        search_model.add_argument(
            '-pdir',
            '--pdir',
            help=('Input directory with models'),
            type=str,
            default=None,
        dest='pdir',
        required=False)
	    
        map_res = parser.add_argument_group()
        map_res.add_argument(
	    '-r',
            '-res',
            '--map_res',
            help=('Resolution of the map'),
            type=float,
            dest='res',
            default=None)
        map_res.add_argument(
            '-r1',
            '-res1',
            '--map_res1',
            help=('Resolution of the map1'),
            type=float,
            dest='res1',
            default=None)
        map_res.add_argument(
            '-r2',
            '-res2',
            '--map_res2',
            help=('Resolution of the map2'),
            type=float,
            dest='res2',
            default=None)
        
        map_threshold = parser.add_argument_group()
        map_threshold.add_argument(
            '-t',
            '--thr',
            help='Map threshold value',
            type=float,
            dest='thr',
            default=None)
        map_threshold.add_argument(
            '-t1',
            '--thr1',
            help='Map1 threshold value',
            type=float,
            dest='thr1',
            default=None)
        map_threshold.add_argument(
            '-t2',
            '--thr2',
            help='Map2 threshold value',
            type=float,
            dest='thr2',
            default=None)
        map_threshold.add_argument(
            '-lT',
            '--tlow',
            help='Lower map threshold value',
            type=float,
            default=None)
        map_threshold.add_argument(
            '-hT',
            '--thigh',
            help='Higher map threshold value',
            type=float,
            default=None)
        
        num_hits = parser.add_argument_group()
        num_hits.add_argument(
            '-hits',
            '--nhits',
            help=('Number of solutions required (1-10000)'),
            type=int,
            dest='num_hits',
            default=100)
        list_inp = parser.add_argument_group()
        list_inp.add_argument(
            '-l',
            '--lfile',
            help=('Input list file'),
            type=str,
            dest='lfile',
            default=100)
        rigid_inp = parser.add_argument_group()
        rigid_inp.add_argument(
            '-rf',
            '--rigidfile',
            help=('Input rigid body file'),
            type=str,
            dest='rigidfile',
            default=None)
        rigid_inp.add_argument(
            '-w',
            '--window',
            help=('Window length'),
            type=int,
            dest='window',
            default=None)
        rigid_out = parser.add_argument_group()
        rigid_out.add_argument(
            '-rigidout',
            '--rigidout',
            help=('Input rigid body file'),
            #type=str,
            dest='rigidout',
            action='store_true',
            default=False)
        
        mode = parser.add_argument_group()
        mode.add_argument(
            '-ref',
            '--mode',
            help=('Input mode'),
            type=int,
            dest='mode',
            default=None)
        apix = parser.add_argument_group()
        apix.add_argument(
            '-s',
            '--apix',
            help=('Input grid spacing'),
            type=float,
            dest='apix',
            default=None)
        
        flag_filt = parser.add_argument_group()
        flag_filt.add_argument(
            '-sf',
            '--sigfac',
            help=('Input Sigma factor for blurring'),
            type=float,
            dest='sigfac',
            default=None)
        flag_filt.add_argument(
            '-nofilt',
            '--nofilt',
            help=('Do not lowpass filter map?'),
            #type=bool,
            dest='nofilt',
            action='store_true',
            default=False)
        flag_filt.add_argument(
            '-bpfilt',
            '--bpfilt',
            help=('bandpass filter map?'),
            #type=bool,
            dest='bpfilt',
            action='store_true',
            default=False)
        flag_filt.add_argument(
            '-mask',
            '--mask',
            help=('Apply mask'),
            #type=bool,
            dest='mask',
            action='store_true',
            default=False)
        flag_filt.add_argument(
            '-smask',
            '--softmask',
            help=('Apply soft mask'),
            #type=bool,
            dest='softmask',
            action='store_true',
            default=False)
        flag_filt.add_argument(
            '-nodust',
            '--nodust',
            help=('Disable dust filter'),
            #type=bool,
            dest='nodust',
            action='store_true',
            default=bool(False))
        flag_filt.add_argument(
            '-softdust',
            '--softdust',
            help=('Disable dust filter'),
            #type=bool,
            dest='softdust',
            action='store_true',
            default=bool(False))
        flag_filt.add_argument(
            '-dp',
            '--dustprob',
            help=('Probability of finding dust among all density particles (greater than this value)'),
            type=float,
            dest='dustprob',
            default=0.2)
        
        flag_shell = parser.add_argument_group()
        flag_shell.add_argument(
            '-noscale',
            '--noscale',
            help=('Disable amplitude scaling'),
            #type=bool,
            dest='noscale',
            action='store_true',
            default=False)
        flag_shell.add_argument(
            '-refscale',
            '--refscale',
            help=('Scale amplitudes based on reference map'),
            #type=bool,
            dest='refscale',
            action='store_true',
            default=False)
        flag_shell.add_argument(
            '-sw',
            '--shellwidth',
            help=('Width of frequency (1/resolution) shell'),
            type=float,
            dest='shellwidth',
            default=0.02)
        
        
        
    
        #
        self.args = parser.parse_args()
        if not self.args.inp_map1 is None and self.args.inp_map2 is None:
            self.args.inp_map = self.args.inp_map1
        elif self.args.inp_map1 is None and not self.args.inp_map2 is None:
            self.args.inp_map = self.args.inp_map2
        if not self.args.pdb1 is None and self.args.pdb2 is None:
            self.args.pdb = self.args.pdb1
        elif self.args.pdb1 is None and not self.args.pdb2 is None:
            self.args.pdb = self.args.pdb2
        if not self.args.thr1 is None and self.args.thr2 is None:
            self.args.pdb = self.args.pdb1
        elif self.args.thr1 is None and not self.args.thr2 is None:
            self.args.thr = self.args.thr2
        '''
            parser.error('-f argument is required in "download" mode.')
        '''

