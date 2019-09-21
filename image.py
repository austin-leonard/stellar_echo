""" Various functions to perform analysis of flare star images from LCO."""
import os
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import cv2
from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from photutils import DAOStarFinder
from PIL import Image
from scipy import signal
from pandas.plotting import autocorrelation_plot


class ImageAnalyze:
    """Analyzes FITS images from the BANZAI pipeline. Should be the reduced (e91) images."""

    def __init__(self, hdul=None, image_data=None, star=None, directory=None,
                 file_list=None, reduced_images=None, n_frames=None, norm=None,
                 positions=None, apertures=None, pixels=None, tar=None, target_center=None,
                 centers=None, list_of_centers=None, orig=None):
        self.hdul = hdul
        self.image_data = image_data
        self.star = star
        self.directory = directory
        self.file_list = file_list
        self.reduced_images = reduced_images
        self.n_frames = n_frames
        self.norm = norm
        self.positions = positions
        self.apertures = apertures
        self.pixels = pixels
        self.tar = tar
        self.target_center = target_center
        self.centers = centers
        self.list_of_centers = list_of_centers
        self.orig = orig

    def load(self):
        """Loads an individual image and a list of reduced fits images to work with."""

        # Load up an image:
        image_file = get_pkg_data_filename(
            "C:/Users/Austin/Desktop/Nanohmics/fits_images/v1396_cyg/ogg0m406-kb27-20190612-0160-e91.fits")
        self.hdul = fits.open(image_file)
        self.image_data = fits.getdata(image_file, ext=0)  # I'm unsure if this line is necessary
        self.star = self.hdul[0].header['object']  # For plot titles and WCS search later

        # Load a directory, and return only the reduced fits images (maybe change in the futures):
        self.directory = os.chdir(
            "C:/Users/Austin/Desktop/Nanohmics/fits_images/v1396_cyg")  # Directory containing fits images
        self.file_list = os.listdir("C:/Users/Austin/Desktop/Nanohmics/fits_images/v1396_cyg")
        self.reduced_images = []  # Contains reduced image data
        self.n_frames = len(self.file_list)

        for j, filename in enumerate(self.file_list):
            if filename.endswith("e91.fits"):
                images = fits.getdata(filename, ext=0)
                self.reduced_images.append(images)
                im = Image.fromarray(self.reduced_images[j]).convert('RGB')
                # im.save('C:/Users/Austin/Desktop/Nanohmics/fits_images/test %s.png' % j) # to save png images

        print(type(self.reduced_images[0]))
        return self.reduced_images

    def find_sources(self, image=0):
        """ Locates the stars in an image using aperture photometry.
            :param image: the image number in the set. 0 is the first image.
            :returns:
            centers: a list of the image coordinates of the center of every star in the image
            apertures: circles surrounding each star
        """
        if image == 0:
            image = self.reduced_images[0]
        elif image != 0:
            index = range(len(self.reduced_images))
            image = self.reduced_images[index[image]]  # this works please don't get rid of this

        data = image.data  # Get the image data. Needs to be an HDUL object
        mean, median, std = sigma_clipped_stats(data,
                                                sigma=3.0)  # Median & sigma of the specified data range.

        daofind = DAOStarFinder(fwhm=3.0, threshold=6. * std)  # Define star criteria
        sources = daofind(data - median)

        self.positions = (sources['xcentroid'], sources['ycentroid'])
        self.apertures = CircularAperture(self.positions, r=4.)
        coords = np.asarray(self.positions)
        xcoords = coords[0]
        ycoords = coords[1]
        self.centers = np.column_stack((xcoords, ycoords))  # Re-format the coordinates in (x,y) fashion
        print(self.centers)
        return self.centers, self.apertures

    def correlate(self, image=0):
        # Performs 2-dimensional correlation with scipy to find stars.
        if image == 0:
            image = self.reduced_images[0]
        elif image != 0:
            index = range(len(self.reduced_images))
            image = self.reduced_images[index]

        # First, define the template star. We'll use the brightest star in the image.

        # Load an image
        img = cv2.imread("C:/Users/Austin/Desktop/Nanohmics/fits_images/test outputs/test 0.png")

        # Convert to gray
        gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

        # Apply a Gaussian blur to the image to reduce noise. Eliminate stray bright pixels.
        gray = cv2.GaussianBlur(gray, (21, 21), 0)
        cv2.h
        # Now call the MinMaxLoc function to find brightest pixel
        (minVal, maxVal, minLoc, maxLoc) = cv2.minMaxLoc(gray)
        print(minVal, maxVal, minLoc, maxLoc)

        # Circle brightest region:
        cv2.circle(gray, maxLoc, 21, (255, 0, 0), 2)

        # So the template star is located at the max location.
        # Now, we'll use the template star as the correlation coefficient.

        # Save the template star location
        max_x = maxLoc[0]
        max_y = maxLoc[1]
        target_image = self.reduced_images[0]
        brightest_pixel = target_image[max_y][max_x]  # Get the location of the brightest pixel in the image

        # Define template star as a range around the brightest pixel
        template = target_image[max_y - 15:max_y + 15, max_x - 15: max_x + 15]

        corr = signal.correlate(target_image, template)

        plt.imshow(corr, norm=self.norm, origin='lower')
        plt.colorbar()
        plt.show()

        # autocorrelation_plot(corr[30])


    def plot_single_image(self, image=0, mark_sources=False, savefig=False):
        """Plots a single image.
                :param image: defaults to the first image in the set.
                :param mark_sources: If True, performs aperture photometry to specify all the sources in the image.
                :param savefig: If True, saves the figure.
                :return: none.
                """
        if image == 0:
            image = self.reduced_images[0]
        elif image != 0:
            index = range(len(self.reduced_images))
            image = self.reduced_images[index[image]]  # THIS FUCKIN WORKS HOLY SHIT

        self.norm = ImageNormalize(vmin=-300, vmax=22000, stretch=SqrtStretch())
        data = image.data

        if mark_sources:
            self.find_sources()

            # self.w_c_s = wcs.WCS(self.hdul[0].header)
            # plt.subplot(projection=w_c_s)

            # Display image with the sources marked
            im = plt.imshow(data, cmap='gray', norm=self.norm, origin='lower', animated=True)
            self.apertures.plot(color='red', lw=1.5, alpha=0.8)
            plt.title("%s, All sources marked" % self.star)
            plt.colorbar()

            if savefig or mark_sources:
                plt.savefig("C:/Users/Austin/Desktop/Nanohmics/fits_images/stars_marked.png", bbox_inches='tight',
                            dpi=1000)
            elif savefig and not mark_sources:
                plt.savefig("C:/Users/Austin/Desktop/Nanohmics/fits_images/image.png", bbox_inches='tight')

            plt.show()
            return im

        elif not mark_sources:
            plt.figure()
            img = plt.imshow(image.data, cmap='gray', norm=self.norm, origin='lower')
            plt.title(f"{self.star} Image")
            plt.colorbar()
            plt.show()
            return img

    def plot_frames(self, save=False):
        """ Loop through all the images in a directory and plot them.
        * NEEDS WORK -- NOT DONE *
        """

        image_directory = self.file_list
        self.norm = ImageNormalize(vmin=-300, vmax=22000, stretch=SqrtStretch())

        fig = plt.figure()  # Get a figure up there
        for ip in image_directory:
            plt.plot(i)
        plt.show()
        plt.xlim(1500, 1550)
        plt.ylim(970, 1020)
        ax = fig.add_subplot(111)  # yeah
        ax.set_title("yep")

        im = ax.imshow(self.reduced_images[0], cmap='gray', norm=self.norm, origin='lower')
        fig.show()
        im.axes.figure.canvas.draw()

        for files in image_directory:
            dat = fits.getdata(files, ext=0)
            ax.set_title(str(files))
            im.set_data(dat)
            if save:
                plt.savefig("C:/Users/Austin/Desktop/Nanohmics/fits_images/uh.png", bbox_inches='tight')
            im.axes.figure.canvas.draw()

    def find_target_center(self, image=0):
        """Locates the center of the target star.
            If no image is specified, defaults to first image in the set.
        """

        if image == 0:
            image = self.reduced_images[0]
        elif image != 0:
            index = range(len(self.reduced_images))
            image = self.reduced_images[index[image]]  # THIS FUCKIN WORKS HOLY SHIT

        # Prepare the image data for WCS transformation
        # image = image.data
        # Going from image coordinates to sky coordinates
        w = wcs.WCS(self.hdul[0].header)
        world = w.wcs_pix2world(self.centers, 0)

        # Get RA and DEC of target star:
        target = SkyCoord.from_name(self.star,
                                    parse=True)  # ICRS RA/DEC, which is what SIMBAD uses. Change to star of interest
        # print("Target star coordinates:\n", "RA:", target.ra.hms, "\n", "DEC:", target.dec)

        # If the RA and DEC of a source in the image sufficiently matches the RA and DEC of the source star from SIMBAD,
        # Print it to the screen.
        for values in world:
            c = SkyCoord(ra=values[0] * u.degree, dec=values[1] * u.degree)
            dif_ra = target.ra.degree - c.ra.degree
            dif_dec = target.dec.degree - c.dec.degree
            self.tar = w.wcs_world2pix(target.ra.degree, target.dec.degree, 0)
            if abs(dif_ra) <= 0.009 and abs(dif_dec) <= 0.009:
                # print("Image coordinates of source:\n", "RA:", c.ra.hms, "within:", abs(dif_ra), "of target", "\n",
                # "DEC:", c.dec, "within:", abs(dif_dec), "of target")
                self.target_center = w.wcs_world2pix(c.ra.degree, c.dec.degree, 0)

                # now grab the pixel coordinates of the image source:
                print("(x,y) pixel location of source center:", self.target_center)
        return self.target_center

    def plot_target_star(self, image=0, save=False):
        """ A plot of just the target star in the image. """

        if image == 0:
            image = self.reduced_images[0]
        elif image != 0:
            index = range(len(self.reduced_images))
            image = self.reduced_images[index[image]]

        img = image.data
        circle = patches.Circle(self.tar, 5, color='green',
                                alpha=0.5)  # Represents real coordinates of the target star.
        # Will be slightly off due to some rounding errors in the calculations.
        self.find_sources()
        fig = plt.figure()
        ax = fig.add_subplot()
        plt.imshow(img, cmap='gray', norm=self.norm)
        plt.title("%s Target Star" % self.star)
        ax.add_patch(circle)
        plt.xlim(self.target_center[0] - 30, self.target_center[0] + 30)
        plt.ylim(self.target_center[1] - 30, self.target_center[1] + 30)
        self.apertures.plot(color='red', lw=1.5, alpha=0.8)
        plt.colorbar()
        if save:
            plt.savefig('C:/Users/Austin/Desktop/Nanohmics/fits_images/target.png', bbox_inches='tight',
                        dpi=1000)
            plt.show()
        else:
            plt.show()

    def get_target_coordinates(self, image=None, directory=None, save=False):
        """Obtain the range of coordinates for the target star in a FITS Image.
        As of right now, only works when plot_single_image has mark_sources=True."""
        global ordered_pair
        if image is None:
            image = self.reduced_images[0]
        elif image is not None:
            image = self.reduced_images[-1]
        if directory is None:
            directory = self.file_list

        self.find_target_center()
        self.plot_target_star()

        # Print the desired pixel range -- the aperture center, plus/minus a few pixels in all directions.

        pix_range_x = [self.target_center[0] - 5, self.target_center[0] + 5]
        pix_range_y = [self.target_center[1] - 5, self.target_center[1] + 5]

        print("Target Star pixel range:\n ", "x:", pix_range_x, "\n", "y:", pix_range_y)
        print("--------------------------------------")
        ordered_pair = np.array((pix_range_x, pix_range_y))
        return ordered_pair, self.pixels

    def track_sway(self):
        """ Track the telescope sway over time by tracking the change in position of one star. """

        # Access the first frame to get the starting (x,y) position.
        print("Starting position: ", self.target_center)

        # Get the location of the target star's center in every image
        self.list_of_centers = []
        for x, image_files in enumerate(self.reduced_images):
            self.find_sources(image=x)
            more_centers = self.find_target_center(image=x)
            self.list_of_centers.append(more_centers)

        # Isolate the x and y values for plotting
        x_values = []
        for a, elements in enumerate(self.list_of_centers):
            new_element = self.list_of_centers[a]
            x = new_element[0]
            x_values.append(x)
        print(x_values)

        y_values = []
        for b, element in enumerate(self.list_of_centers):
            n_element = self.list_of_centers[b]
            y = n_element[1]
            y_values.append(y)
        print(y_values)

        # Make a plot
        first_element = self.list_of_centers[0]
        first_x = first_element[0]
        first_y = first_element[1]

        last_element = self.list_of_centers[-1]
        last_x = last_element[0]
        last_y = last_element[-1]

        new_x = x_values[1:8]
        new_y = y_values[1:8]

        plt.scatter(first_x, first_y, color='red')
        plt.scatter(last_x, last_y, color='green')
        plt.scatter(new_x, new_y, color='blue')
        plt.xlim(1500, 1530)
        plt.ylim(980, 1010)
        plt.title("Change in position of center over time")
        plt.show()


if __name__ == '__main__':
    i = ImageAnalyze()
    i.load()
    i.correlate()
    # i.plot_single_image(image=0, mark_sources=False)
    # i.find_sources(image=2)
    # i.find_target_center(image=2)
    # i.plot_target_star(image=2)
    # i.plot_frames()
    # i.get_target_coordinates()
    # i.track_sway()
    # i.plot_centers_over_time()

'''
Things to add:
- Find a way to implement 2d correlation to find stars more effectively. 
- Implement 2x2 boxcar binning of images 
- Pixel masks and stuff
'''
