function load_array(arrayPtr::Ptr, length::Integer, juliaArray)
    Threads.@threads for i = 1:length
        a= unsafe_load(arrayPtr, i)
        juliaArray[i] = isnan(a) ? 0 : a
    end
    print("finished loading after ", length, " elements\n")
end

function store_array(arrayPtr::Ptr, length::Integer, juliaArray)
    juliaArray = reshape(juliaArray, (:,1)) 
    if length <= size(juliaArray)[1]
        Threads.@threads for i = 1:length
            unsafe_store!(arrayPtr, juliaArray[i], i)
        end
    else
        print("Array lengths do not match", size(juliaArray)[1], " vs ", length, "\n")
    end
    print("finished storing after ", length, " elements, make sure this is the intended number of stored entries\n")
end

function store_array(arrayPtr::Ptr, juliaArray::AbstractArray, startHeight::Integer=1)
    size = Base.size(juliaArray)
    start = 0
    if Base.length(size) == 2
        length = size[1]*size[2]
        start = (startHeight-1)*size[1]+1
    elseif Base.length(size) == 3
        length = size[1]*size[2]*size[3]
        start = (startHeight-1)*size[1]*size[2]+1
        print("Starting at ", start, "\n")
        print("Loading ", length, " elements\n")
    elseif Base.length(size) == 4
        length = size[1]*size[2]*size[3]*size[4]
        start = (startHeight-1)*size[1]*size[2]*size[3]+1
        print("Starting at ", start, "\n")
        print("Loading ", length, " elements\n")
    else
        error("Invalid size of array")
    end
    juliaArray = reshape(juliaArray, (:,1))
    Threads.@threads for i = 1:length
        unsafe_store!(arrayPtr, juliaArray[i], start+i-1)
    end
    print("finished storing after ", length, " elements, make sure this is the intended number of stored entries\n")
    return
end

