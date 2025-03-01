# Utilities to predict using enformer
# Author: Temi
# Date: Thurs Feb 2 2023

import functools

@functools.lru_cache(5)
def get_model(model_path, dl_package='tensorflow'):
    """
    Return a tensorflow model

    Parameters:
        model_path: str
            A path to where the tensorflow model exists
    Returns: 
        a tensorflow model
    """

    if dl_package == 'tensorflow':
        import tensorflow as tf
        return tf.saved_model.load(model_path).model
    elif dl_package == 'pytorch':
        import torch
        import enformer_pytorch
        device = torch.device("cpu")

        stored_state = torch.load(model_path, map_location=device)

        model = enformer_pytorch.Enformer()

        # model.load_state_dict(stored_state["model_state_dict"])
        model.load_state_dict(stored_state)

        model.to('cuda')
  
        return model

def enformer_predict_on_sequence(model, sample_input, dl_package='tensorflow'):
    """
    given a compatible sequence that has been one-hot encoded, predict on ENFORMER

    Parameters:
        model: a tensorflow model
        sample_input: a (1, 393216, 4) np.array that is a one-hot encoding of a sequence

    Returns: A dictionary
        of the form {'haplotype': _predictions_}
        _predictions_ is a numpy array of shape (17, 5313) numpy array of predictions
    """
    
    prediction_output = {}
    for haplotype, sequence_encoding in sample_input.items():
        if not sequence_encoding.shape == (1, 393216, 4):
            raise Exception(f'[ERROR] Fatal. Input sequence shape is not appropriate')
        # prediction = model.predict_on_batch(sequence_encoding)['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]

        if dl_package == 'tensorflow':
            import tensorflow as tf
            prediction = model.predict_on_batch(sequence_encoding)['human'].numpy()
        elif dl_package == 'pytorch':
            import torch
            prediction_raw = model(torch.tensor(sequence_encoding).to("cuda"))
            prediction = prediction_raw['human'].cpu().detach().numpy()
            del prediction_raw
            torch.cuda.empty_cache()
            print("PREDICTION COMPLETE")

        prediction_output[haplotype] = prediction
        
    return(prediction_output)


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]