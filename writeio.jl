module WriteIO

export WriteFile

function WriteFile(headpart::Function, bodypart::Function, tailpart::Function,
                   filename="test.txt")
  if isfile(filename)
    println("file ", filename, " already, exist, replacing it ...")
    rm(filename)
  end
  filepointer = open(filename, "a")
  pr(x="") = println(filepointer, x)

  headpart(pr);
  bodypart(pr);
  tailpart(pr);

  println("Done generating ", filename)
  close(filepointer)
end

end
