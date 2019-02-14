import DelayEmbeddings: Dataset

"""
    cembed(data::Dataset)

Returns an embedding consisting of a zero-lagged, unmodified version of `data`.
"""
function cembed(data::Dataset)
    if size(data, 1) > size(data, 2)
        #info("Treating each row as a point")
    	dim = size(data, 2)
        which_pos = [i for i = 1:dim]
        which_lags = [0 for i in 1:dim]
        cembed([data[:, i] for i = 1:dim], which_pos, which_lags)
    else
    	#info("Treating each column of data as a point")
    	dim = size(data, 1)
        which_pos = [i for i = 1:dim]
        which_lags = [0 for i in 1:dim]
        cembed([data[i, :] for i = 1:dim], which_pos, which_lags)
    end
end

"""
    cembed(d::Dataset,
        which_pos::Vector{Int},
        which_lags::Vector{Int})

Returns a state space embedding of the columns of `d`.

## Arguments
* `d::Dataset`: The columns of `d` contains the data series to use
    for the embedding.
* `which_pos::Vector{Int}``. The length of in_which_pos gives the dimension
    of the embedding. The value of the ith element of `in_which_pos`
    indicates which column of `d` goes in the i-th column of the embedding.
* `which_lags::Vector{Int}` sets the lag in each column of the reconstrution.
    Must be the same length as `dimension(d)`
    - **Example**: if `in_which_pos = [2, 2, 1]`, then
        `at_what_lags = [1, 0, -1]` means that the lag in column 1 is 1, the
        lag in the second column is 0 and the lag in the third column is -1.
"""
function cembed(data::Dataset,
        which_pos,
        which_lags)

    if size(data, 1) > size(data, 2)
        #info("Treating each row as a point")
    	dim = size(data, 2)
    	cembed([data[:, i] for i = 1:dim], which_pos, which_lags)
    else
    	#info("Treating each column of data as a point")
    	dim = size(data, 1)
        cembed([data[i, :] for i = 1:dim], which_pos, which_lags)
    end
end

function Embedding(data::Dataset,
        which_pos,
        which_lags)
    if size(data, 1) > size(data, 2)
        #info("Treating each row as a point")
        dim = size(data, 2)
        cembed([data[:, i] for i = 1:dim], which_pos, what_lags)
    else
        #info("Treating each column of data as a point")
        dim = size(data, 1)
        cembed([data[i, :] for i = 1:dim], which_pos, what_lags)
    end
end

export cembed