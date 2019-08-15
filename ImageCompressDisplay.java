import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;

import static java.lang.Math.*;
import java.util.Arrays;

public class ImageCompressDisplay {
    int width = 512;
	int height = 512;

	JFrame frame = new JFrame();

    JLabel lbIm1 = new JLabel();
	JLabel lbIm2 = new JLabel();
    
    JLabel Text1 = new JLabel();
    JLabel Text2 = new JLabel();

    GridBagLayout GLayout = new GridBagLayout();

	BufferedImage OriginalImage;
    BufferedImage DCTImage;
    BufferedImage DWTImage;

    static double[][] cosineMatrix = new double[8][8];

    //Orignal of R, G, B
    int[][] RMatrix = new int[height][width];
    int[][] GMatrix = new int[height][width];
    int[][] BMatrix = new int[height][width];

    //DCT of R, G, B
    int[][] DCT_RMatrix = new int[height][width];
    int[][] DCT_GMatrix = new int[height][width];
    int[][] DCT_BMatrix = new int[height][width];
        
    //IDCT of R, G, B
    int[][] IDCT_RMatrix = new int[height][width];
    int[][] IDCT_GMatrix = new int[height][width];
    int[][] IDCT_BMatrix = new int[height][width];
        
    //DWT of R, G, B
    double[][] DWT_RMatrix = new double[height][width];
    double[][] DWT_GMatrix = new double[height][width];
    double[][] DWT_BMatrix = new double[height][width];

    //IDWT of R, G, B
    int[][] IDWT_RMatrix = new int[height][width];
    int[][] IDWT_GMatrix = new int[height][width];
    int[][] IDWT_BMatrix = new int[height][width];

	/** Read Image RGB
	 *  Reads the image of given width and height at the given imgPath into the provided BufferedImage.
	 */
	private void readImageRGB(int width, int height, String imgPath, BufferedImage img)
	{
		try
		{
        
			int frameLength = width*height*3;

			File file = new File(imgPath);
			RandomAccessFile raf = new RandomAccessFile(file, "r");
			raf.seek(0);

			long len = frameLength;
			byte[] bytes = new byte[(int) len];

			raf.read(bytes);
             
			int ind = 0;
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width; x++)
				{
					int r = bytes[ind];
					int g = bytes[ind+height*width];
					int b = bytes[ind+height*width*2]; 

                    //Turn sign int to unsign int.
                    RMatrix[y][x] = r & 0xff;
                    GMatrix[y][x] = g & 0xff;
                    BMatrix[y][x] = b & 0xff;
                
                    /*
                    RMatrix[y][x] = r;
                    GMatrix[y][x] = g;
                    BMatrix[y][x] = b;
                    */

					int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
					//int pix = ((a << 24) + (r << 16) + (g << 8) + b);
					img.setRGB(x,y,pix);
					ind++;
				}
			}
		}
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	public void showIms(String[] args){

		// Read in the specified image
        String ImagePath = args[0];

        //Initial BufferedImage.
        OriginalImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        DCTImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        DWTImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        //Read Source Image
        readImageRGB(width, height, ImagePath, OriginalImage);
        
        int ncoef = Integer.parseInt(args[1]);
        
        //Initiate consine transform values in matrix
        for(int i=0; i<8; i++){
            for(int j=0; j<8; j++){
                cosineMatrix[i][j] = cos((2*i+1)*j*3.1415926/16);
            }
        }
        
        if(ncoef!=-1){
            int mcoef = ncoef/4096;
            
            //Discret Cosine Transform
            DCTTransform(mcoef);

            //Do Inverse Discret Cosine Transform
            IDCTTransform();

            DWTTransform(ncoef);

            IDWTTransform();

            showDCTDWT(0);
        }
        else{
            int iter = 1;
            int current = 4096;
            while(current <= 262144){
                ncoef = current;
                int mcoef = ncoef/4096;

                DCTTransform(mcoef);

                IDCTTransform();

                DWTTransform(ncoef);

                IDWTTransform();
            
                //showDCTDWT(0);

                try{
                    Thread.sleep(1000);
                }
                catch(InterruptedException ex){
                    Thread.currentThread().interrupt();
                }               

                showDCTDWT(iter);
                iter++;

                if(current == 262144){
                    current = 0;
                    iter = 1;
                }               
                current+=4094;
            }
        }


		// Use label to display the image
		/*
        frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		lbIm1 = new JLabel(new ImageIcon(OriginalImage));

		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 0;
		c.gridy = 0;

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 1;
		frame.getContentPane().add(lbIm1, c);

		frame.pack();
		frame.setVisible(true);
        */
	}

    private void showDCTDWT(int iteration){
	    //frame = new JFrame();

        //lbIm1 = new JLabel();
	    //lbIm2 = new JLabel();
    
        //Text1 = new JLabel();
        //Text2 = new JLabel();

        //GLayout = new GridBagLayout();
        
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                int r = IDCT_RMatrix[i][j];
                int g = IDCT_GMatrix[i][j];
                int b = IDCT_BMatrix[i][j];

                int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
				//int pix = ((a << 24) + (r << 16) + (g << 8) + b);
	            DCTImage.setRGB(j, i, pix);

                int ir = (int)IDWT_RMatrix[i][j];
                int ig = (int)IDWT_GMatrix[i][j];
                int ib = (int)IDWT_BMatrix[i][j];

                int ipix = 0xff000000 | ((ir & 0xff) << 16) | ((ig & 0xff) << 8) | (ib & 0xff);
				//int pix = ((a << 24) + (r << 16) + (g << 8) + b);
	            DWTImage.setRGB(j, i, ipix);
            }
        }

        frame.getContentPane().setLayout(GLayout);
        if(iteration != 0){
            Text1.setText("Discret Cosine Transform Iteration " + iteration);
            Text2.setText("Discret Wavelet Transform Iteration " + iteration);
        }   
        else{
            Text1.setText("Discret Cosine Transform");
            Text2.setText("Discret Wavelet Transform");
        }
        
        Text1.setHorizontalAlignment(SwingConstants.CENTER);
        Text2.setHorizontalAlignment(SwingConstants.CENTER);

        lbIm1.setIcon(new ImageIcon(DCTImage));
        lbIm2.setIcon(new ImageIcon(DWTImage));

    	GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 0;
		c.gridy = 0;
        frame.getContentPane().add(Text1, c);

        c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 1;
		c.gridy = 0;
        frame.getContentPane().add(Text2, c);

    	c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.gridx = 0;
		c.gridy = 1;
        frame.getContentPane().add(lbIm1, c);

    	c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.gridx = 1;
		c.gridy = 1;
        frame.getContentPane().add(lbIm2, c);

		frame.pack();
		frame.setVisible(true);
    }

    public double[][] zigzag(double[][] matrix, int mcoef){
        //int length = matrix.length-1;
        int length = matrix.length;
        int count = 1;            

        double[][] reference = new double[length*2-1][length];
        int[][][] reference_location = new int[length*2-1][length][2];
        int[] counter = new int[length*2-1];
        
        for(int i=1; i<=(2*length)-1; i++){
            int start_col = Math.max(0, i-length);
            int record = Math.min(i, length - start_col);
            record = Math.min(record, length);
            for(int j=0; j<record; j++){
                reference[i-1][counter[i-1]] = matrix[Math.min(length, i)-j-1][start_col+j];
                reference_location[i-1][counter[i-1]][0] = Math.min(length, i) - j - 1;
                reference_location[i-1][counter[i-1]][1] = start_col+j; 
                counter[i-1] = counter[i-1] + 1;
            }
        }
        
        int flag = 0;
        for(int i=0; i<length*2-1; i++){
            if(flag%2 == 0){
                for(int j=0; j<counter[i]; j++){
                    if(count > mcoef){
                        matrix[reference_location[i][j][0]][reference_location[i][j][1]] = 0;
                        count++;
                    }
                    else{
                        count++;
                    }
                }
            }
            else{
                for(int j=counter[i]-1; j>=0; j--){
                     if(count > mcoef){
                        matrix[reference_location[i][j][0]][reference_location[i][j][1]] = 0;
                        count++;
                    }
                    else{
                        count++;
                    }
                }
            }
            flag++;
        }
 
        return matrix;
    }

    private void DWTDecomposition(double[][] DWT_Matrix, int[][] Matrix, int ncoef){
        //System.out.println("height: " + height);
        //System.out.println("width: " + width);

        double[][] Buffer_Matrix = new double[height][width];
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                Buffer_Matrix[i][j] = Matrix[i][j];
            }
        }

        //System.out.println("Buffer_Matrix: " + Buffer_Matrix.length);
        //System.out.println("Buffer_Matrix[0]: " + Buffer_Matrix[0].length);

        for(int row=0; row<width; row++){
            Buffer_Matrix[row] = Decomposition(Buffer_Matrix[row]);
        }
        
        //System.out.println("Buffer_Matrix: " + Buffer_Matrix.length);
        //System.out.println("Buffer_Matrix[0]: " + Buffer_Matrix[0].length);
        
        Buffer_Matrix = Transpose(Buffer_Matrix);
        
        for(int col=0; col<height; col++){
            Buffer_Matrix[col] = Decomposition(Buffer_Matrix[col]);
        }

        Buffer_Matrix = Transpose(Buffer_Matrix);
        
        Buffer_Matrix = zigzag(Buffer_Matrix, ncoef);

        
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                DWT_Matrix[i][j] = Buffer_Matrix[i][j];
            }
        }
        
        //return Buffer_Matrix;
        //DWT_Matrix = Buffer_Matrix;
    }

    private double[] Decomposition(double[] List){
        int length = List.length;
        while(length > 0){
            //List = DecompositionProcess(List, length);
            double[] Buffer_List = Arrays.copyOf(List, List.length);
            for(int i=0; i<length/2; i++){
                Buffer_List[i] = (List[i*2] + List[i*2+1])/2;
                Buffer_List[length/2 + i] = (List[2*i] - List[2*i+1])/2;
            } 
            List = Arrays.copyOf(Buffer_List, Buffer_List.length);

            length/=2;
        
        }
        return List;
    }

    private double[][] Transpose(double[][] Matrix){
        double[][] Buffer = new double[height][width];
        //System.out.println(Buffer.length);
        //System.out.println(Buffer[0].length);
        //System.out.println(Matrix.length);
        //System.out.println(Matrix[0].length);

        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                Buffer[i][j] = Matrix[j][i];
            }
        }
        return Buffer;
    }

    private void IDWTComposition(int[][] IDWT_Matrix, double[][] DWT_Matrix){
        int[][] Buffer_Matrix = new int[height][width];

        DWT_Matrix = Transpose(DWT_Matrix);
        for(int col=0; col<height; col++){
            DWT_Matrix[col] = Composition(DWT_Matrix[col]);
        }
        DWT_Matrix = Transpose(DWT_Matrix);

        for(int row=0; row<width; row++){
            DWT_Matrix[row] = Composition(DWT_Matrix[row]);
        }

        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                Buffer_Matrix[i][j] = (int)Math.round(DWT_Matrix[i][j]);
                if(Buffer_Matrix[i][j] < 0){
                    Buffer_Matrix[i][j] = 0;
                }
                else if(Buffer_Matrix[i][j] > 255){
                    Buffer_Matrix[i][j] = 255;
                }
            }
        }
        //return Buffer_Matrix;
        
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                IDWT_Matrix[i][j] = Buffer_Matrix[i][j];
            }
        } 
        
        //DWT_Matrix = (int)Buffer_Matrix;
    }

    private double[] Composition(double[] List){
        int index = 1;
        while(index <= List.length){
            //List = CompositionProcess(List, index);
            double[] Buffer_List = Arrays.copyOf(List, List.length);
         
            for(int i=0; i<index/2; i++){
                Buffer_List[2*i] = List[i] + List[index/2 + i];
                Buffer_List[2*i + 1] = List[i] - List[index/2 + i];
            }    
            List = Arrays.copyOf(Buffer_List, Buffer_List.length);

            index*=2;
        }
        return List;
    }

    private void DWTTransform(int ncoef){
        //Three chanel R, G, B
        DWTDecomposition(DWT_RMatrix, RMatrix, ncoef);
        DWTDecomposition(DWT_GMatrix, GMatrix, ncoef);
        DWTDecomposition(DWT_BMatrix, BMatrix, ncoef);
    }

    private void IDWTTransform(){
        //Three chanel R, G, B
        IDWTComposition(IDWT_RMatrix, DWT_RMatrix);
        IDWTComposition(IDWT_GMatrix, DWT_GMatrix);
        IDWTComposition(IDWT_BMatrix, DWT_BMatrix);
    }

    private void DCTTransform(int mcoef){
         //Store the block values
        double[][] RBlock = new double[8][8];
        double[][] GBlock = new double[8][8];
        double[][] BBlock = new double[8][8];

        for(int i=0; i<height; i+=8){
            for(int j=0; j<width; j+=8){
                for(int k=0; k<8; k++){
                    for(int l=0; l<8; l++){
                        double cu = 1;
                        double cv = 1;
                        double RResult = 0;
                        double GResult = 0;
                        double BResult = 0;
                            
                        if(k == 0){
                            cu = 0.707;
                        }
                        if(l == 0){
                            cv = 0.707;
                        }
                        
                        for(int m=0; m<8; m++){
                            for(int n=0; n<8; n++){
                                int R = (int)RMatrix[i+m][j+n];
                                int G = (int)GMatrix[i+m][j+n];
                                int B = (int)BMatrix[i+m][j+n];
                                RResult += R*cosineMatrix[m][k]*cosineMatrix[n][l];
                                GResult += G*cosineMatrix[m][k]*cosineMatrix[n][l];
                                BResult += B*cosineMatrix[m][k]*cosineMatrix[n][l];

                            }
                        }
                        
                        RBlock[k][l] = (int)Math.round(RResult*0.25*cu*cv);
                        GBlock[k][l] = (int)Math.round(GResult*0.25*cu*cv);
                        BBlock[k][l] = (int)Math.round(BResult*0.25*cu*cv);
                    }
                }
                
                RBlock = zigzag(RBlock, mcoef);
                GBlock = zigzag(GBlock, mcoef);
                BBlock = zigzag(BBlock, mcoef);
                
                for(int k=0; k<8; k++){
                    for(int l=0; l<8; l++){
                        DCT_RMatrix[i+k][j+l] = (int)RBlock[k][l];
                        DCT_GMatrix[i+k][j+l] = (int)GBlock[k][l];
                        DCT_BMatrix[i+k][j+l] = (int)BBlock[k][l];
                    }
                }
            }
        }
    }


    public void IDCTTransform(){
        for(int i=0; i<height; i+=8){
            for(int j=0; j<width; j+=8){
                for(int k=0; k<8; k++){
                    for(int l=0; l<8; l++){
                        double dR = 0, dG = 0, dB = 0;
                        
                        for(int m = 0; m<8; m++){
                            for(int n=0; n<8; n++){
                                double cu = 1, cv = 1;
                                
                                double R, G, B;
                                if(m == 0){
                                    cu = 0.707;
                                    
                                }
                                if(n == 0){
                                    cv = 0.707;
                                }
                                
                                R = DCT_RMatrix[i+m][j+n];
                                G = DCT_GMatrix[i+m][j+n];
                                B = DCT_BMatrix[i+m][j+n];
                                
                                dR += cu*cv*R*cosineMatrix[k][m]*cosineMatrix[l][n];
                                dG += cu*cv*G*cosineMatrix[k][m]*cosineMatrix[l][n];
                                dB += cu*cv*B*cosineMatrix[k][m]*cosineMatrix[l][n];
                                
                            }
                        }

                        dR = dR*0.25;
                        dG = dG*0.25;
                        dB = dB*0.25;
                        //Filter out the out of range data.
                        //If it is over 255 to be 255.
                        //If it is under 0 to be 0.
                        if(dR <= 0){
                            dR = 0;
                        }
                        else if(dR >= 255){
                            dR = 255;
                        }
                        if(dG <= 0){
                            dG = 0;
                        }
                        else if(dG >= 255){
                            dG = 255;
                        } 
                        if(dB <= 0){
                            dB = 0;
                        }
                        else if(dB >= 255){
                            dB = 255;
                        }

                        IDCT_RMatrix[i+k][j+l] = (int)dR;
                        IDCT_GMatrix[i+k][j+l] = (int)dG;
                        IDCT_BMatrix[i+k][j+l] = (int)dB;
                    }
                }
            }
        }
    }

	public static void main(String[] args) {
		ImageCompressDisplay ren = new ImageCompressDisplay();
		ren.showIms(args);
	}

}
