import cv2
import numpy as np
import os
import matplotlib.pyplot as plt

inputs = "C:/Users/martk/Documents/wur/masterthesis/microscopy/multiphoton/tiff"
threshold = 80 # Can be changed when not as much change is observed in sample(s).
files = os.listdir(inputs)
files = sorted([file for file in files if file.endswith(".tif")])
with open(f'Ablation_analysis_second_threshold2.csv', 'w') as f:

    if not os.path.exists('Output'):
        os.mkdir('Output')

    for file in files:
        print(file)

        folder = f'Output/{file[:-4]}'
        if not os.path.exists(folder):
            os.mkdir(folder)
        images = cv2.imreadmulti(os.path.join(inputs, file), [], cv2.IMREAD_UNCHANGED)[1]
       # images = cv2.imreadmulti(file, [], cv2.IMREAD_UNCHANGED)[1]
        images = [image.astype(np.int16) for image in images]
        maxlist = []
        for index, image in enumerate(images[:-1]):
            imagesDiff = images[index+1] - images[index]
            imagesDiff = cv2.GaussianBlur(imagesDiff, (3, 3), 0)
            maxlist += [np.min(imagesDiff),np.max(imagesDiff)]
            images[index] = imagesDiff
        print("after for-loop")
        #SumProjection = np.sum(np.abs(images)[:-1], axis=0)
        SumProjection = np.mean(np.abs(images)[:-1], axis=0)
        #SumProjection = np.max(np.abs(images)[:-1], axis=0)
        #cv2.imwrite(f'{folder}_SumProjection.tif',SumProjection)
        cv2.imwrite(f'{folder}_SumProjection.tif', SumProjection.astype(np.uint16))

        mask = np.zeros_like(SumProjection).astype(np.uint8)
        mask[SumProjection > threshold] = 255
        cv2.imwrite(f'{folder}_Mask.tiff', mask)

        contours, hierarchy = cv2.findContours(mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        contours = sorted(contours, key=cv2.contourArea, reverse=True)
        contours = [contours[i] for i in range(len(contours)) if cv2.contourArea(contours[i]) > 1000]

        ContourMask = np.zeros_like(mask)
        cv2.drawContours(ContourMask, contours, 0, 255, -1)
        cv2.imwrite(f'{folder}_Mask_contour.tiff', ContourMask)
        absmax = max(np.abs([min(maxlist), max(maxlist)]))
        f.write(str(file[:-4])+ ',')
        for index, image in enumerate(images[:-1]):
            #DO CALCULATIONS HERE, WITH np.abs(image[mask == 255])
            value = np.mean(np.abs(image[ContourMask == 255]))
            f.write(str(value) + ',')

            pltMask = np.zeros_like(image)
            pltMask[ContourMask == 255] = image[ContourMask == 255]
            plt.matshow(pltMask, cmap=plt.cm.bwr,vmin=-absmax, vmax=absmax)
            plt.colorbar()
            plt.savefig(f"{folder}/{index}.png")
            plt.close()
        f.write('\n')
f.close()
