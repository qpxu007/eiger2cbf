/*
EIGER HDF5 to CBF converter
 Written by Takanori Nakane

To build:

Linux

gcc -std=c99 -o eiger2cbf -g \
 -I$HOME/prog/dials/modules/cbflib/include \
 -L$HOME/prog/dials/build/lib -Ilz4 \
 eiger2cbf.c \
 lz4/lz4.c lz4/h5zlz4.c \
 bitshuffle/bshuf_h5filter.c \
 bitshuffle/bshuf_h5plugin.c \
 bitshuffle/bitshuffle.c \
 /usr/lib/x86_64-linux-gnu/libhdf5_hl.a \
 /usr/lib/x86_64-linux-gnu/libhdf5.a \
 -lcbf -lm -lpthread -lz -ldl

Mac OS X

gcc -std=c99 -o eiger2cbf -g \
 -ICBFlib-0.9.5.2/include -Ilz4 \
 eiger2cbf.c \
 lz4/lz4.c lz4/h5zlz4.c \
 bitshuffle/bshuf_h5filter.c \
 bitshuffle/bshuf_h5plugin.c \
 bitshuffle/bitshuffle.c \
 $HOME/local/lib/libhdf5_hl.a \
 $HOME/local/lib/libhdf5.a \
 CBFlib-0.9.5.2/lib/libcbf.a \
 -lm -lpthread -lz -ldl

*/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "getopt.h"
#include "unistd.h"
#include "omp.h"

#include "cbf.h"
#include "cbf_simple.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "omp.h"

extern const H5Z_class2_t H5Z_LZ4;
extern const H5Z_class2_t bshuf_H5Filter;
void register_filters()
{
  H5Zregister(&H5Z_LZ4);
  H5Zregister(&bshuf_H5Filter);
}

// Function to extract the filename from a full path
const char *extractFilename(const char *path)
{
  const char *filename = strrchr(path, '/'); // Find the last occurrence of '/'

  if (filename == NULL)
  {
    // If no '/' found, return the full path
    return path;
  }
  else
  {
    // Return the filename (skip the '/')
    return filename + 1;
  }
}

int main(int argc, char **argv)
{
  int xpixels = -1, ypixels = -1, beamx = -1, beamy = -1, nimages = -1, depth = -1, countrate_cutoff = -1, ntrigger = 1;
  int from = -1, to = -1;
  double pixelsize = -1, wavelength = -1, distance = -1, count_time = -1, frame_time = -1, osc_width = -1, osc_start = -9999, thickness = -1;
  char detector_sn[256] = {}, description[256] = {}, version[256] = {};
  char renumber = 1;

  hid_t hdf;

  fprintf(stderr, "EIGER HDF5 to CBF converter (version 160530 gmca mod)\n");
  fprintf(stderr, " written by Takanori Nakane\n");
  fprintf(stderr, " see https://github.com/biochem-fan/eiger2cbf for details.\n\n");

  int opt;
  char *prefix = NULL;
  while ((opt = getopt(argc, argv, "s:e:p:xh")) != -1)
  {
    switch (opt)
    {
    case 's':
      from = atoi(optarg);
      break;
    case 'e':
      to = atoi(optarg);
      break;
    case 'p':
      prefix = optarg;
      break;
    case 'x':
      renumber = -1; // disable renumbering
      fprintf(stderr, "renumbering based on angle disabled\n");
      break;
    case 'h':
      fprintf(stderr, "Usage: %s -s start -e end -p prefix master_file\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  char *master_file = argv[optind];
  if (master_file == NULL || access(master_file, F_OK) == -1)
  {
    fprintf(stderr, "Usage: %s -s start -e end -p prefix master_file\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  printf("master file: %s\n", master_file);

  if (prefix == NULL)
  {
    const char *filename = extractFilename(master_file);
    char *location = strstr(filename, "master.");

    int prefix_len = location - filename;

    prefix = (char *)realloc(prefix, (prefix_len + 1) * sizeof(char));
    // char *prefix[prefix_len + 1];
    strncpy(prefix, filename, prefix_len);
    prefix[prefix_len] = '\0';
  }
  printf("Prefix: %s\n", prefix);

  register_filters();

  hdf = H5Fopen(master_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (hdf < 0)
  {
    fprintf(stderr, "Failed to open file %s\n", master_file);
    return -1;
  }

  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/nimages", &nimages);
  herr_t status = H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/ntrigger", &ntrigger);

  if (from == -1 && to == -1)
  {
    printf("nimages=%d ntrigger=%d total_images=%d\n", nimages, ntrigger, ntrigger * nimages);
    from = 1;
    to = nimages * ntrigger;
  }
  else if (from >= 1 && to == -1)
  {
    to = from;
  }

  nimages = nimages * ntrigger;

  fprintf(stderr, "Going to convert frame %d to %d.\n", from, to);

  H5Eset_auto(0, NULL, NULL); // Comment out this line for debugging.

  fprintf(stderr, "Metadata in HDF5:\n");
  H5LTread_dataset_string(hdf, "/entry/instrument/detector/description", description);
  fprintf(stderr, " /entry/instrument/detector/description = %s\n", description);
  H5LTread_dataset_string(hdf, "/entry/instrument/detector/detector_number", detector_sn);
  fprintf(stderr, " /entry/instrument/detector/detector_number = %s\n", detector_sn);
  H5LTread_dataset_string(hdf, "/entry/instrument/detector/detectorSpecific/software_version", version);
  fprintf(stderr, " /entry/instrument/detector/detectorSpecific/software_version = %s\n", version);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/bit_depth_image", &depth);
  if (depth > 0)
  {
    fprintf(stderr, " /entry/instrument/detector/bit_depth_image = %d\n", depth);
  }
  else
  {
    fprintf(stderr, " WARNING: /entry/instrument/detector/bit_depth_image is not avaialble. We assume 16 bit.\n");
    depth = 16;
  }
  unsigned int error_val = (unsigned int)(((unsigned long long)1 << depth) - 1);

  // Saturation value
  // Firmware >= 1.5
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/saturation_value", &countrate_cutoff);
  if (countrate_cutoff > 0)
  {
    fprintf(stderr, " /entry/instrument/detector/detectorSpecific/saturation_value = %d\n", countrate_cutoff);
  }
  else
  {
    // Firmware >= 1.4
    fprintf(stderr, "  /entry/instrument/detector/detectorSpecific/saturation_value not present. Trying another place.\n");
    H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff", &countrate_cutoff);
    if (countrate_cutoff > 0)
    {
      fprintf(stderr, " /entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff = %d\n", countrate_cutoff);
      countrate_cutoff++;
    }
    else
    {
      fprintf(stderr, "  /entry/instrument/detector/detectorSpecific/countrate_correction_count_cutoff not present. Trying another place.\n");
      // < 1.4
      H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/detectorModule_000/countrate_correction_count_cutoff", &countrate_cutoff);
      if (countrate_cutoff > 0)
      {
        fprintf(stderr, " /entry/instrument/detector/detectorSpecific/detectorModule_000/countrate_correction_count_cutoff = %d\n", countrate_cutoff);
        fprintf(stderr, "  WARNING: The use of this field is not recommended now.\n");
        fprintf(stderr, "  You might want to change the OVERLOAD setting in your subsequent processing.\n");
        countrate_cutoff++;
      }
      else
      {
        fprintf(stderr, "  /entry/instrument/detector/detectorSpecific/detectorModule_000/countrate_correction_count_cutoff not present.\n");
        countrate_cutoff = error_val - 1;
        fprintf(stderr, "  As a last resort, we will put an arbitrary large number (%d) in the header.\n", countrate_cutoff);
        fprintf(stderr, "  You might want to change the OVERLOAD setting in your subsequent processing.\n");
      }
    }
  }

  H5LTread_dataset_double(hdf, "/entry/instrument/detector/sensor_thickness", &thickness); // in m
  if (thickness > 0)
  {
    fprintf(stderr, " /entry/instrument/detector/sensor_thickness = %f (um)\n", thickness * 1E6);
  }
  else
  {
    thickness = 450E-6;
    fprintf(stderr, " /entry/instrument/detector/sensor_thickness is not avaialble. We assume it is %f um\n", thickness * 1E6);
  }
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/x_pixels_in_detector", &xpixels);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/y_pixels_in_detector", &ypixels);
  fprintf(stderr, " /entry/instrument/detector/detectorSpecific/{x,y}_pixels_in_detector = (%d, %d) (px)\n",
          xpixels, ypixels);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_x", &beamx);
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_y", &beamy);
  fprintf(stderr, " /entry/instrument/detector/beam_center_{x,y} = (%d, %d) (px)\n", beamx, beamy);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/count_time", &count_time); // in m
  fprintf(stderr, " /entry/instrument/detector/count_time = %f (sec)\n", count_time);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/frame_time", &frame_time); // in
  fprintf(stderr, " /entry/instrument/detector/frame_time = %f (sec)\n", frame_time);
  H5LTread_dataset_double(hdf, "/entry/instrument/detector/x_pixel_size", &pixelsize); // in
  fprintf(stderr, " /entry/instrument/detector/x_pixel_size = %f (m)\n", pixelsize);

  // Detector distance

  H5LTread_dataset_double(hdf, "/entry/instrument/detector/distance", &distance); // Firmware >= 1.7
  if (distance > 0)
  {
    fprintf(stderr, " /entry/instrument/detector/distance = %f (m)\n", distance);
  }
  else
  {
    fprintf(stderr, "  /entry/instrument/detector/distance not present. Trying another place.\n");

    H5LTread_dataset_double(hdf, "/entry/instrument/detector/detector_distance", &distance); // Firmware< 1.7
    if (distance > 0)
    {
      fprintf(stderr, " /entry/instrument/detector/detector_distance = %f (m)\n", distance);
    }
    else
    {
      fprintf(stderr, "  /entry/instrument/detector/detector_distance not present.\n");
      fprintf(stderr, " WARNING: detector distance was not defined! \"Detector distance\" field in the output is set to -1.\n");
    }
  }

  // Wavelength
  H5LTread_dataset_double(hdf, "/entry/sample/beam/incident_wavelength", &wavelength); // Firmware >= 1.7
  if (wavelength > 0)
  {
    fprintf(stderr, " /entry/sample/beam/incident_wavelength = %f (A)\n", wavelength);
  }
  else
  {
    fprintf(stderr, "  /entry/sample/beam/incident_wavelength not present. Trying another place.\n");

    H5LTread_dataset_double(hdf, "/entry/instrument/beam/wavelength", &wavelength);
    if (wavelength > 0)
    {
      fprintf(stderr, " /entry/instrument/beam/wavelength = %f (A)\n", wavelength);
    }
    else
    {
      fprintf(stderr, "  /entry/instrument/beam/wavelength not present. Trying another place.\n");

      H5LTread_dataset_double(hdf, "/entry/instrument/monochromator/wavelength", &wavelength);
      if (wavelength > 0)
      {
        fprintf(stderr, " /entry/instrument/monochromator/wavelength = %f (A)\n", wavelength);
      }
      else
      {
        fprintf(stderr, "  /entry/instrument/monochromator/wavelength not present. Trying another place.\n");

        H5LTread_dataset_double(hdf, "/entry/instrument/beam/incident_wavelength", &wavelength); // Firmware 1.6
        if (wavelength > 0)
        {
          fprintf(stderr, " /entry/instrument/beam/incident_wavelength = %f (A)\n", wavelength);
        }
        else
        {
          fprintf(stderr, "  /entry/instrument/beam/incident_wavelength not present.\n");
        }
      }
    }
  }
  if (wavelength < 0)
  {
    fprintf(stderr, " WARNING: wavelength was not defined! \"Wavelength\" field in the output is set to -1.\n");
  }

  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega_range_average", &osc_width);
  if (osc_width > 0)
  {
    fprintf(stderr, " /entry/sample/goniometer/omega_range_average = %f (deg)\n", osc_width);
  }
  else
  {
    fprintf(stderr, " WARNING: oscillation width was not defined. \"Start_angle\" field in the output is set to 0!\n");
    osc_width = 0;
  }
  unsigned int *buf = (unsigned int *)malloc(sizeof(unsigned int) * xpixels * ypixels);
  signed int *buf_signed = (signed int *)malloc(sizeof(signed int) * xpixels * ypixels);
  signed int *pixel_mask = (signed int *)malloc(sizeof(signed int) * xpixels * ypixels);
  if (buf == NULL || buf_signed == NULL || pixel_mask == NULL)
  {
    fprintf(stderr, "Failed to allocate image buffer.\n");
    return -1;
  }

  // TODO: Is it always in omega?
  int bufsize = (nimages < 100000) ? 100000 : nimages;
  double *angles = (double *)malloc(bufsize * sizeof(double));
  // I don't know why but nimages can be too small ...
  if (angles == NULL)
  {
    fprintf(stderr, "failed to allocate buffer for omega.\n");
    return -1;
  }
  angles[0] = -9999;
  H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega", angles);
  fprintf(stderr, "\n");

  hid_t entry, group;
  entry = H5Gopen2(hdf, "/entry", H5P_DEFAULT);
  if (entry < 0)
  {
    fprintf(stderr, "/entry does not exist!\n");
    return -1;
  }

  pixel_mask[0] = -9999;
  H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/pixel_mask", pixel_mask);
  if (pixel_mask[0] == -9999)
  {
    fprintf(stderr, "WARNING: failed to read the pixel mask from /entry/instrument/detector/detectorSpecific/pixel_mask.\n");
    fprintf(stderr, " Thus, we mask pixels whose intensity is %u (= (2 ^ bit_depth_image) - 1) by converting them to -1. \n", error_val);
    fprintf(stderr, " However, this might mask overloaded (saturated) pixels as well.\n");
  }

  // Check if /entry/data present
  group = H5Gopen2(entry, "data", H5P_DEFAULT);
  if (group < 0)
  {
    group = entry; // leak!
  }

  int block_start = 1;
  if (H5LTfind_dataset(group, "data_000000"))
  {
    fprintf(stderr, "This dataset starts from data_000000.\n");
    block_start = 0;
  }
  else
  {
    fprintf(stderr, "This dataset starts from data_000001.\n");
  }

  char data_name[20] = {};
  hid_t data, dataspace;
  int number_per_block = 0;

  // Open the first data block to get the number of frames in a block
  snprintf(data_name, 20, "data_%06d", block_start);
  data = H5Dopen2(group, data_name, H5P_DEFAULT);
  dataspace = H5Dget_space(data);
  if (data < 0)
  {
    fprintf(stderr, "failed to open /entry/%s\n", data_name);
    return -1;
  }
  if (H5Sget_simple_extent_ndims(dataspace) != 3)
  {
    fprintf(stderr, "Dimension of /entry/%s is not 3!\n", data_name);
    return -1;
  }
  hsize_t dims[3];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  number_per_block = dims[0];
  fprintf(stderr, "The number of images per data block is %d.\n", number_per_block);

  H5Sclose(dataspace);
  H5Dclose(data);

  fprintf(stderr, "\nFile analysis completed.\n\n");
  int frame;

#pragma omp parallel for num_threads(32) private(hdf, data, dataspace, data_name, frame, osc_start) shared(angles, buf_signed, buf)
  for (frame = from; frame <= to; frame++)
  {
    fprintf(stderr, "Converting frame %d (%d / %d)\n", frame, frame - from + 1, to - from + 1);
    if (angles[0] != -9999)
    {
      osc_start = angles[frame - 1];
      fprintf(stderr, " /entry/sample/goniometer/omega[%d] = %.3f %.2f(1-indexed)\n", frame, osc_start, angles[frame-1]);
    }
    else
    {
      fprintf(stderr, " oscillation start not defined. \"Start_angle\" field in the output is set to 0!\n");
      osc_start = osc_width * frame; // old firmware
    }

    int modified_frame = frame;
    if (renumber == 1 && osc_width >= 1e-6)
    {
      modified_frame = (int)round((osc_start - angles[0]) / osc_width + 1);
    }

    char filename[4096];
    snprintf(filename, 4096, "%s%06d.cbf", prefix, modified_frame);
    fprintf(stderr, "frame=%i --> %i  osc=%.3f outfile=%s\n", frame, modified_frame, osc_start, filename);

    if (frame > nimages)
    {
      fprintf(stderr, "WARNING: invalid frame number specified. %d is bigger than nimages (%d)\n", frame, nimages);
      // Due to a firmware bug, nimages can be smaller than the actual value.
      // So we don't exit here
    }

    char header_format[] =
        "\n"
        "# Detector: %s, S/N %s\n"
        "# Pixel_size %de-6 m x %de-6 m\n"
        "# Silicon sensor, thickness %.6f m\n"
        "# Exposure_time %f s\n"
        "# Exposure_period %f s\n"
        "# Count_cutoff %d counts\n"
        "# Wavelength %f A\n"
        "# Detector_distance %f m\n"
        "# Beam_xy (%d, %d) pixels\n"
        "# Start_angle %f deg.\n"
        "# Angle_increment %f deg.\n";

    char header_content[4096] = {};
    snprintf(header_content, 4096, header_format,
             description, detector_sn,
             thickness,
             (int)(pixelsize * 1E6), (int)(pixelsize * 1E6),
             count_time, frame_time, countrate_cutoff, wavelength, distance,
             beamx, beamy, osc_start, osc_width);

    // Now open the required data

    int block_number = block_start + (frame - 1) / number_per_block;
    int frame_in_block = (frame - 1) % number_per_block;
    //    fprintf(stderr, " frame %d is in data_%06d frame %d (1-indexed).\n",
    //            frame, block_number, frame_in_block + 1);

    snprintf(data_name, 20, "data_%06d", block_number);
    char err_msg[4096] = "";
    data = H5Dopen2(group, data_name, H5P_DEFAULT);
    dataspace = H5Dget_space(data);
    if (data < 0)
    {
      sprintf(err_msg, "failed to open /entry/%s\n", data_name);
    }
    if (H5Sget_simple_extent_ndims(dataspace) != 3)
    {
      sprintf(err_msg, "Dimension of /entry/%s is not 3!\n", data_name);
    }

    // Get the frame
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    hsize_t offset_in[3] = {frame_in_block, 0, 0};
    hsize_t offset_out[3] = {0, 0, 0};
    hsize_t count[3] = {1, ypixels, xpixels};
    hid_t memspace = H5Screate_simple(3, dims, NULL);
    if (memspace < 0)
    {
      sprintf(err_msg, "failed to create memspace\n");
    }

    int ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset_in, NULL,
                                  count, NULL);
    if (ret < 0)
    {
      sprintf(err_msg, "select_hyperslab for file failed\n");
    }
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL,
                              count, NULL);
    if (ret < 0)
    {
      sprintf(err_msg, "select_hyperslab for memory failed\n");
    }
    ret = H5Dread(data, H5T_NATIVE_UINT, memspace, dataspace, H5P_DEFAULT, buf);
    if (ret < 0)
    {
      sprintf(err_msg, "H5Dread for image failed. Wrong frame number? frame=%d\n", frame);
    }

    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Dclose(data);

    if (strlen(err_msg) > 0)
    {
      fprintf(stderr, "--Error--: %s", err_msg);
    }

    /////////////////////////////////////////////////////////////////
    // Reading done. Here output starts...

    FILE *fh = fopen(filename, "wb");

    // create a CBF
    cbf_handle cbf;
    cbf_make_handle(&cbf);
    cbf_new_datablock(cbf, "image_1");

    // put a miniCBF header
    cbf_new_category(cbf, "array_data");
    cbf_new_column(cbf, "header_convention");
    cbf_set_value(cbf, "SLS_1.0");
    cbf_new_column(cbf, "header_contents");
    cbf_set_value(cbf, header_content);

    // put the image
    cbf_new_category(cbf, "array_data");
    cbf_new_column(cbf, "data");
    int i;
    for (i = 0; i < xpixels * ypixels; i++)
    {
      if ((pixel_mask[0] != -9999 && pixel_mask[i] == 1) || // the pixel mask is available
          (pixel_mask[0] == -9999 && buf[i] == error_val))
      { // not available
        buf_signed[i] = -1;
      }
      else if (pixel_mask[0] != -9999 && pixel_mask[i] > 1)
      { // the pixel mask is 2, 4, 8, 16
        buf_signed[i] = -2;
      }
      else
      {
        buf_signed[i] = buf[i];
      }
    }
    cbf_set_integerarray_wdims_fs(cbf,
                                  CBF_BYTE_OFFSET,
                                  1, // binary id
                                  buf_signed,
                                  sizeof(int),
                                  1, // signed?
                                  xpixels * ypixels,
                                  "little_endian",
                                  xpixels,
                                  ypixels,
                                  0,
                                  0); // padding
    cbf_write_file(cbf, fh, 1, CBF, MSG_DIGEST | MIME_HEADERS | PAD_4K, 0);
    // no need to fclose() here as the 3rd argument "readable" is 1
    cbf_free_handle(cbf);
  }

  H5Gclose(group);
  H5Fclose(hdf);

  free(buf);
  free(buf_signed);
  free(angles);

  fprintf(stderr, "\nAll done!\n");

  return 0;
}
