import cv2
import numpy as np
import os

inputs = "C:/Users/martk/Documents/wur/masterthesis/phenotyping/TIFFfilesACT-ADF-phenotyping"
files = os.listdir(inputs)
files = [file for file in files if file.endswith(".tif")]

for file in files:
    print(file)
    image = cv2.imread(os.path.join(inputs,file))
    mask = np.zeros_like(image[:, :, 0])
    mask[image[:, :, 0] < 80] = 255

    # cv2.imwrite(f'{file[:-4]}_test0.tiff', mask)

    contours, hierarchy = cv2.findContours(mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)

    cutoff = 9
    minsize = cv2.contourArea(contours[cutoff])
    contours = [contours[i] for i in range(len(contours)) if cv2.contourArea(contours[i]) > 70]

    mask2 = np.zeros_like(image)
    for i in range(len(contours[:cutoff])):
        cv2.drawContours(mask2, contours, i, [255,0,0], -1)  # Draw filled contour in mask
        area = cv2.contourArea(contours[i])
        #print(f'intitial area = {area}')

        for index2, contour2 in enumerate(contours):
            if i != index2:
                m = cv2.moments(contour2)
                try:
                    cX = int(m["m10"] / m["m00"])
                    cY = int(m["m01"] / m["m00"])
                    if cv2.pointPolygonTest(contours[i], (cX, cY), False) == -1:
                        if cv2.contourArea(contours[index2]) < minsize:  # and cv2.contourArea(contours[index2]) > 100:
                            cv2.drawContours(mask2, contours, index2, [0, 0, 0], -1)  # Draw filled contour in mask
                            area = area - cv2.contourArea(contours[index2])
                except:
                    ""
        print(f'{area}')

    cv2.imwrite(f'{file[:-4]}_test.tiff', mask2)
