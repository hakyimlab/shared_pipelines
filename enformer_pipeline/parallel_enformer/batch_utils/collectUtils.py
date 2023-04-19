import numpy as np

def parse_bins_and_tracks(unparsed_bins, unparsed_tracks):
    def split_range(range_str):
        if "-" in range_str:
            spl_range_str = [int(x) for x in range_str.split("-")]
            return [x for x in range(spl_range_str[0],spl_range_str[1]+1)]
        else:
            return [int(range_str)]
    
    if unparsed_bins == -1:
        bins_indices = None
    else:
        bins_indices = []
        bins_split = unparsed_bins.split(",")
        for b in bins_split:
            bins_indices += split_range(b)
        bins_indices.sort()
    
    if unparsed_tracks == -1:
        tracks_indices = None
    else:
        tracks_indices = []
        tracks_split = unparsed_tracks.split(",")
        for t in tracks_split:
            tracks_indices += split_range(t)
        tracks_indices.sort()

        # with open(mapping_path,"w") as mp:
        #     mp.write("stored_track_index\tmodel_track_index\n")
        #     for i,x in enumerate(tracks_indices):
        #         mp.write("\t".join([str(i),str(x)])+"\n")

    return bins_indices,tracks_indices

def collect_bins_and_tracks(predictions, bins_indices, tracks_indices):
    
    if bins_indices == None:
        if tracks_indices == None:
            return(predictions[:, :])
        else:
            return(predictions[:, tracks_indices])
    else:
        if tracks_indices == None:
            return(predictions[bins_indices, :])
        else:
            return(predictions[bins_indices, tracks_indices])
        


# import numpy as np
# x = np.array(((1,2,3),(4,5,6)))

# b_inds = "0-1"
# t_inds = "2-2"

# bins_indices,tracks_indices = parse_bins_and_tracks(b_inds,t_inds)
# out = collect_bins_and_tracks(x,bins_indices,tracks_indices)
# print(out)