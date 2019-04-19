import os
import logging

logger = logging.getLogger(__name__)


def cleanup_files(animation_func):
    def inner(*args,**kwargs):
        animation_func(*args,**kwargs)
        files = args[1]
        for f in files:
            os.remove(f)
    return inner

class Animation(object):
    """
    Adopted from https://zulko.wordpress.com/2012/09/29/animate-your-3d-plots-with-pythons-matplotlib/

    """

    def __init__(self,movie_name = 'movie.gif', movie_path = None):
        self.anim_filename = movie_name
        self.anim_format = movie_name.split('.')[1]
        self.anim_path = movie_path if movie_path \
                else os.path.join(os.getcwd(),self.anim_filename)


    @staticmethod
    def make_3Dviews(ax,angles,elevation=0, width=4, height = 8,
                prefix='tmprot_',**kwargs):
        """
        Makes jpeg pictures of the given 3d ax, with different angles.
        Args:
            ax (3D axis): te ax
            angles (list): the list of angles (in degree) under which to
                           take the picture.
            width,height (float): size, in inches, of the output images.
            prefix (str): prefix for the files created.

        Returns: the list of files created (for later removal)
        """

        files = []
        ax.figure.set_size_inches(width,height)

        for i,angle in enumerate(angles):

            ax.view_init(elev = elevation, azim=angle)
            fname = '%s%03d.jpeg'%(prefix,i)
            ax.figure.savefig(fname)
            files.append(fname)

        return files

    # Transform a series of picture into animation
    @cleanup_files
    def make_movie(self,files,fps=10,bitrate=1800,**kwargs):
        """
        Uses mencoder, produces a .mp4/.ogv/... movie from a list of
        picture files.
        """

        output_name, output_ext = os.path.splitext(self.animpath)
        command = { '.mp4' : 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc lavc\
                             -lavcopts vcodec=msmpeg4v2:vbitrate=%d'
                             %(",".join(files),fps,output_name,bitrate)}

        command['.ogv'] = command['.mp4'] + \
            '; ffmpeg -i %s.mp4 -r %d %s'%(output_name,fps,self.anim_path)

        logger.debug(command[output_ext])
        output_ext = os.path.splitext(self.anim_path)[1]
        os.system(command[output_ext])

#    @cleanup_files
    def make_gif(self,files,delay=100, repeat=True,**kwargs):
        """
        Uses imageMagick to produce an animated .gif from a list of
        picture files.
        """

        loop = -1 if repeat else 0
        os.system('convert -delay %d -loop %d %s %s'
                  %(delay,loop," ".join(files),self.anim_path))

    @cleanup_files
    def make_strip(self,files,**kwargs):
        """
        Uses imageMagick to produce a .jpeg strip from a list of
        picture files.
        """
        os.system('montage -tile 1x -geometry +0+0 %s %s'\
                  %(" ".join(files),self.anim_path))



