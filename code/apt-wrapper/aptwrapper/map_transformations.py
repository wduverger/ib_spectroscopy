from typing import Tuple, Union
import numpy as np
from skimage.registration import phase_cross_correlation
import xarray as xr
from scipy import interpolate

from scipy.ndimage import shift
from skimage.transform import ProjectiveTransform, SimilarityTransform, rotate, warp, EuclideanTransform
from skimage.measure import ransac
from skimage.feature import plot_matches, ORB, match_descriptors, plot_matches

from .plot_maps import pixel_size_xy, force_2d_coords, plot_map
import matplotlib.pyplot as plt

## Matrix interpolation
#######################

def get_3d_idx(hmap):
    i = hmap.ypix.data
    j = hmap.xpix.data
    r = np.stack([i, j, np.ones_like(i)], axis=0)
    return r
def get_3d_coords(hmap):
    i = hmap.Y.data
    j = hmap.X.data
    r = np.stack([i, j, np.ones_like(i)], axis=0)
    return r
def get_img_transform(hmap):
    y = hmap.Y.data
    x = hmap.X.data
    m = np.array([[y[1,0] - y[0,0], x[1,0] - x[0,0], y[0,0]],
                  [y[0,1] - y[0,0], x[0,1] - x[0,0], x[0,0]],
                  [0              , 0              , 1     ]])
    return m
def idx_to_coords(idx, m):
    coords = np.einsum('ij,jkl->ikl', m, idx)
    # return coords
    y, x = coords[0], coords[1]
    return (y, x)
def coords_to_idx(coords, m):
    idx = np.einsum('ij,jkl->ikl', np.linalg.inv(m), coords)
    # return idx
    i, j = idx[0], idx[1]
    return (i, j)

def new_interp(hmap, grid_x, grid_y):
    m = get_img_transform(hmap)
    coords = np.stack([grid_y, grid_x, np.ones_like(grid_y)], axis=0)
    idx = coords_to_idx(coords, m)

    x_ = xr.DataArray(idx[1], dims=['y_', 'x_'])
    y_ = xr.DataArray(idx[0], dims=['y_', 'x_'])

    return (
        hmap.interp(x=x_, y=y_)
        .drop_vars(['x', 'y'])  # drop old index values
        .rename(x_='x', y_='y')  # rename temporary dimensions
        # when extrapolating, make sure z values are np.nan, but X and Y are still defined
        .assign_coords(X=(('y', 'x'), grid_x), Y=(('y', 'x'), grid_y))
    )


## Matrix transforms
####################

def rotate(hmap, deg, **kwargs):
    t = np.deg2rad(deg)
    transform = np.array([
        [ np.cos(t), np.sin(t), 0],
        [-np.sin(t), np.cos(t), 0],
        [0,          0,         1],
    ])

    return apply_warp(hmap, transform, **kwargs)

def translate(hmap, dx_um, dy_um, pixels=False, inverted=False):

    if pixels:
        pixel_size_x = np.abs(hmap.X[0,1]-hmap.X[0,0]).item()
        pixel_size_y = np.abs(hmap.Y[1,0]-hmap.Y[0,0]).item()
        dy_um = -pixel_size_y * dy_um
        dx_um =  pixel_size_x * dx_um

    if inverted:
        dy_um = -dy_um
        dx_um = -dx_um


    transform = np.array([
        [1, 0, float(dy_um)],
        [0, 1, float(dx_um)],
        [0, 0, 1           ],
    ])
    # return apply_warp(hmap, transform, **kwargs)
    return (
        hmap
        .assign_coords(X=(('y','x'), hmap.X.data + dx_um))
        .assign_coords(Y=(('y','x'), hmap.Y.data + dy_um))
    )

def apply_warp(hmap, transform, inverted=False, pixels=False, origin='center'):
    """ Apply a matrix transformation to a map."""

    # Copy the warp matrix so we don't modify the original
    transform = transform.copy()

    # Invert the matrix if requested
    if inverted: 
        transform = np.linalg.inv(transform)

    # If pixels==True, then convert the matrix to pixel units
    if pixels:
        pixel_size_x = np.abs(hmap.X.diff('x')[0,0]).item()
        pixel_size_y = np.abs(hmap.Y.diff('y')[0,0]).item()
        transform[0, -1] = -pixel_size_y * transform[0, -1]
        transform[1, -1] =  pixel_size_x * transform[1, -1]

    # Define origin point
    if origin == 'center':
        pivot_x = hmap.X.mean().item()
        pivot_y = hmap.Y.mean().item()
    elif origin == 'upper left':
        pivot_x = hmap.X.min().item()
        pivot_y = hmap.Y.max().item()
    elif origin == 'lower left':
        pivot_x = hmap.X.min().item()
        pivot_y = hmap.Y.min().item()
    else:
        pivot_x, pivot_y = 0, 0

    # Translate the map so that the pivot point is at the origin
    hmap = translate(hmap, -pivot_x, -pivot_y)

    # Construct the coordinate matrix with dimensions compatible with the warp matrix (3xNxN)
    xyt0 = np.stack([hmap.Y, hmap.X, np.ones_like(hmap.Y)], axis=0)

    # Apply the warp matrix
    xyt1 = np.einsum('ij,jkl', transform, xyt0)

    # Extract and save the new coordinates
    hmap = hmap.assign_coords(
        Y=(('y', 'x'), xyt1[0]),
        X=(('y', 'x'), xyt1[1]), 
    )

    return translate(hmap, pivot_x, pivot_y)

def invert_warp(mat):
    mat = mat.copy()
    mat[0:2, 0:2] = np.linalg.inv(mat[0:2, 0:2])
    mat[0:2, 2] *= -1
    return mat


def interp_like(hmap1, hmap2):
    """
    Interpolate hmap1 onto the coordinates of hmap2
    """
    
    hmap1 = force_2d_coords(hmap1)
    hmap2 = force_2d_coords(hmap2)
    return new_interp(hmap1, hmap2.X.to_numpy(), hmap2.Y.to_numpy())


## Simple drift measurement
###########################

def transform_coordinate_space(r: Union[Tuple, np.ndarray]) -> np.ndarray:
    return np.array(r)[::-1] * np.array([1, -1])

def calculate_drift(hmap1, hmap2, mask=None):
    """ 
    Use phase cross-correlation to find how much hmap2 is translated with
    respect to hmap1. Returns an array of x,y drift 
    """
    assert (pixel_size_xy(hmap1) == pixel_size_xy(hmap2)).all()
    assert hmap1.shape == hmap2.shape

    r = phase_cross_correlation(hmap2.to_numpy(), hmap1.to_numpy(), reference_mask=mask)
    drift_pixels = r if mask is not None else r[0]

    # Convert to data coords and um
    drift_um = transform_coordinate_space(drift_pixels) * pixel_size_xy(hmap1)
    return drift_um

def calculate_drift_speed(hmap1, hmap2):
    dx = calculate_drift(hmap1, hmap2)
    dt = hmap2['TimeStamp'] - hmap1['TimeStamp']

    return dx/dt.item()*1e9

## Binary operations
####################

def bin_op(map1, map2, ref1, ref2, fun):
    dx = calculate_drift(ref1, ref2)
    map2 = translate(map2, -dx[0], -dx[1])
    map2 = interp_like(map2, map1)

    return fun(map1, map2)

# Older implementation
# def interp_xarr(hmap: xr.DataArray, xx, yy, method='linear') -> xr.DataArray:

#     desired_shape = xx.shape

#     # Actual
#     xx0 = hmap.X.data
#     yy0 = hmap.Y.data
#     xii0 = np.broadcast_to(hmap.x.data, xx0.shape)
#     yii0 = np.broadcast_to(hmap.y.data, yy0.shape).T

#     xx0, yy0 = xx0.flatten(), yy0.flatten()
#     xii0, yii0 = xii0.flatten(), yii0.flatten()

#     # Interpolated (linear interp slows things down but is necessary if you
#     # don't want jagged edges)
#     xii = interpolate.griddata((xx0, yy0), xii0, (xx.flatten(), yy.flatten()), method=method)
#     yii = interpolate.griddata((xx0, yy0), yii0, (xx.flatten(), yy.flatten()), method=method)

#     xii = xii.reshape(desired_shape)
#     yii = yii.reshape(desired_shape)

#     x_ = xr.DataArray(xii, dims=['y_', 'x_'])
#     y_ = xr.DataArray(yii, dims=['y_', 'x_'])

#     return (
#         hmap.interp(x=x_, y=y_)
#         .drop_vars(['x', 'y'])  # drop old index values
#         .rename(x_='x', y_='y')  # rename temporary dimensions
#         # when extrapolating, make sure z values are np.nan, but X and Y are still defined
#         .assign_coords(X=(('y', 'x'), xx), Y=(('y', 'x'), yy))
#     )


## Image registration
#####################

def apply_skew(hmap, drift_speed):
    x, y = np.meshgrid(hmap.x, hmap.y)
    nx, ny = len(hmap.x), len(hmap.y)
    linetime = 1/hmap['Tags.ScanRate'].item()

    hmap = hmap.assign_coords(
        t=(('y', 'x'), (ny-y)*linetime+(nx-x)*linetime/nx)
    )

    return (
        hmap
        .assign_coords(X=(('y', 'x'), (hmap.X+drift_speed[0]*hmap.t).transpose('y', 'x').data))
        .assign_coords(Y=(('y', 'x'), (hmap.Y+drift_speed[1]*hmap.t).data))
    )

from matplotlib.patches import Rectangle

def inner_bounding_box(hmap: xr.DataArray, plot=False) -> Tuple[float, float, float]:
    """ 
    Returns the center and the size (edge len) of a square that fits inside
    a skewed heightmap. 
    """
    x_min = hmap.X.min(dim='x').max().item()
    x_max = hmap.X.max(dim='x').min().item()
    y_min = hmap.Y.min(dim='y').max().item()
    y_max = hmap.Y.max(dim='y').min().item()

    d = min(y_max-y_min, x_max-x_min)
    cx = (x_min+x_max)/2
    cy = (y_min+y_max)/2

    if plot:
        plot_map(hmap)
        plt.gca().add_patch(Rectangle(
            (cx-d/2, cy-d/2), d, d,
            color='k', fill=False
        ))

    return cx, cy, d


def register_orb(hmap1, hmap2, min_samples=10, residual_threshold=1, max_trials=300):
    orb = ORB(n_keypoints=500, fast_threshold=0.05)
    orb.detect_and_extract(hmap1)
    keypoints1 = orb.keypoints
    descriptors1 = orb.descriptors
    orb.detect_and_extract(hmap2)
    keypoints2 = orb.keypoints
    descriptors2 = orb.descriptors
    matches12 = match_descriptors(descriptors1, descriptors2, cross_check=True)
    src = keypoints2[matches12[:, 1]][:, ::-1]
    dst = keypoints1[matches12[:, 0]][:, ::-1]
    model_robust, _ = ransac(
        (src, dst), 
        EuclideanTransform,
        min_samples=min_samples,
        residual_threshold=residual_threshold,
        max_trials=max_trials,
    )

    if np.isnan(model_robust.params).any():
        raise "ORB returned nan. No transformation found."
    return model_robust

def apply_orb(hmap, model):
    m = np.copy(model.params)
    m[0:2,-1] = m[0:2,-1] * hmap.Y.diff('y').mean().item()
    m[0:2,-1] = [m[1, -1], m[0, -1]]
    m

    return apply_warp(hmap, m)
