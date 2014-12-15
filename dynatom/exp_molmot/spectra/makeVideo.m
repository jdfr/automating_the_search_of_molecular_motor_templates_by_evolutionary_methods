function makeVideo(nameInput, nameVideo)

system(['C:\utils\ffmpeg\bin\ffmpeg -y -i ' nameInput ' -vcodec libx264 -level 41 -crf 23 -g 250 -coder 1 ' nameVideo]);