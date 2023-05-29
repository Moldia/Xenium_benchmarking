import squidpy as sq
import numpy as np
from squidpy.im._container import ImageContainer
from squidpy.im._segment import SegmentationModel
from typing import Union,  Optional, Any, Mapping, Callable, Sequence, TYPE_CHECKING, Tuple
from squidpy._utils import NDArrayA

def segment_nuclei(
    img: ImageContainer,
    layer: Optional[str] = None,
    library_id: Union[str, Sequence[str], None] = None,
    method: Union[str, SegmentationModel, Callable[..., NDArrayA]] = "watershed",
    channel: Optional[int] = 0,
    chunks: Union[str, int, Tuple[int, int], None] = None,
    lazy: bool = False,
    layer_added: Optional[str] = None,
    copy: bool = False,
    **kwargs: Any,
) -> Optional[ImageContainer]:
    """Squidpy segment wrapper function
    Based on https://github.com/scverse/squidpy version 1.2.2
    This function will also smooth the image via ``process``

    Parameters
    ----------
    img : ImageContainer
        High-resolution image.
    layer : Optional[str], optional
        Image layer in `img` that should be processed. If None and only 1 layer is present, 
        it will be selected., by default None
    library_id : Union[str, Sequence[str], None], optional
        Name of the Z-dimension(s) that this function should be applied to. 
        For not specified Z-dimensions, the identity function is applied. 
        If None, all Z-dimensions are segmented separately, by default None
    method : Union[str, SegmentationModel, Callable[..., NDArrayA]], optional
        Segmentation method to use. Valid options are:
        ``watershed`` - skimage.segmentation.watershed().
        Alternatively, any callable() can be passed as long as it has the following signature: 
        ``numpy.ndarray (height, width, channels)`` -> ``numpy.ndarray (height, width[, channels])``,
        by default "watershed"
    channel : Optional[int], optional
        Channel index to use for segmentation. If None, use all channels, by default 0
    chunks : Union[str, int, Tuple[int, int], None], optional
        Number of chunks for dask. For automatic chunking, use ``chunks = 'auto'``, by default None
    lazy : bool, optional
        Whether to lazily compute the result or not. Only used when ``chunks != None``, by default False
    layer_added : Optional[str], optional
        Layer of new image layer to add into img object. If None, use 'segmented_{model}'., by default None
    copy : bool, optional
        If True, return the result, otherwise save it to the image container, by default False

    Returns
    -------
    Optional[ImageContainer]
        If copy = True, returns a new container with the segmented image in '{layer_added}'.
        Otherwise, modifies the img with the following key:
        ``squidpy.im.ImageContainer ['{layer_added}']``
    """  
      

    if (kwargs is not None) and ("blur_std" in kwargs):
        if kwargs["blur_std"] > 0:
            sq.im.process(
                img, 
                layer="image", 
                method="smooth", #skimage.filters.gaussian, 
                layer_added="image",
                sigma=kwargs["blur_std"],
                truncate=4.0,
            )
        del kwargs["blur_std"]
            
    return sq.im.segment(img=img, layer= "image", library_id=library_id, method=method, 
        channel=channel, chunks=chunks, lazy=lazy, layer_added=layer_added, copy=copy, **kwargs)


def segment_cellpose(
    img: NDArrayA, 
    hyperparams: Optional[dict]
) -> NDArrayA:
    """Run cellpose and get masks

    Parameters
    ----------
    img : NDArrayA
        Can be list of 2D/3D images, or array of 2D/3D images, or 4D image array
    Returns
    -------
    NDArray
        labelled image, where 0=no masks; 1,2,...=mask labels
    """
    from cellpose import models
    
    # Set model type
    if (hyperparams is not None) and ("model_type" in hyperparams):
        model_type = hyperparams["model_type"]
    else:
        model_type = 'nuclei'
    
    # Init model
    model = models.Cellpose(model_type=model_type)
    
    # Predict
    if hyperparams is not None:
        if "model_type" in hyperparams:
            del hyperparams["model_type"]
        res, _, _, _ = model.eval(
            img,
            channels=[0, 0],
            **hyperparams
        )
    else:
        res, _, _, _ = model.eval(
            img,
            channels=[0, 0]
        )
        
    return res

def segment_binning(
    img: NDArrayA,
    bin_size: int
) -> NDArrayA:

    # Get shape of image
    n = np.shape(img)[0]
    m = np.shape(img)[1]

    # Create grids of coordinates, and combine to form bins
    x = np.floor(np.mgrid[0:n, 0:m][0] / bin_size)
    y = np.floor(np.mgrid[0:n, 0:m][1] / bin_size)
    bins = x*(np.ceil(m/bin_size)) + y + 1

    return bins
