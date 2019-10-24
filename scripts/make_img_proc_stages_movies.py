import os
import imageutils as iu

def get_identifiers_slicebounds(sector):

    if sector == 6:
        projid = 1500 # initial project id for this sector
    elif sector == 7:
        projid = 1516
    elif sector == 8:
        projid = 1532
    elif sector == 9:
        projid = 1548
    elif sector == 10:
        projid = 1564

    identifiers = []
    slicebounds = []

    for cam in range(1,5):
        for ccd in range(1,5):

            identifiers.append(
                (sector, cam, ccd, projid)
            )

            if projid == 1500: # orion b here
                slicebounds.append(
                    [slice(1,513), slice(300,812)]
                )
            else:
                slicebounds.append(
                    [slice(300,812), slice(300,812)]
                )

            projid += 1

    return identifiers, slicebounds


def main(sector=6):

    identifiers, slicebounds = get_identifiers_slicebounds(sector)

    basedir = '/nfs/phtess2/ar0/TESS/FFI/PROJ/IMG_PROC_STAGES'
    moviedir = '/nfs/phtess2/ar0/TESS/FFI/MOVIES'

    for i,s in zip(identifiers, slicebounds):

        outdir = os.path.join(
            basedir, 'sector{}_cam{}_ccd{}_projid{}'.format(
                i[0], i[1], i[2], i[3]
            )
        )
        if not os.path.exists(outdir):
            print('made {}'.format(outdir))
            os.mkdir(outdir)

        iu.plot_stages_of_img_proc_sector_cam_ccd(sector=i[0], cam=i[1], ccd=i[2],
                                                  projid=i[3], overwrite=0,
                                                  outdir=outdir, slicebounds=s)

        imgglob = os.path.join(outdir, 'tess2*_img_proc_stages.png')

        outmp4path = os.path.join(
            moviedir, 'img_proc_stages_sector{}_cam{}_ccd{}_projid{}.mp4'.format(
                i[0], i[1], i[2], i[3]
            )
        )
        if not os.path.exists(outmp4path):
            iu.make_mp4_from_jpegs(imgglob, outmp4path,
                                   ffmpegpath='/home/lbouma/bin/ffmpeg',
                                   verbose=True)

        else:
            print('found {} and skipped'.format(outmp4path))


if __name__ == "__main__":

    sector = 6

    main(sector=sector)
