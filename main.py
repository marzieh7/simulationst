import argparse
from preprocess_tumor import preprocess_tumor
from preprocess_visium_spatial import preprocess_visium_spatial

import os
import torch
def main(args):

    # Device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # Load and preprocess data
    
    # مرحله ۱: پیش‌پردازش
    #preprocess_tumor(args.sc_tumor, args.preprocessed_sc)

    preprocess_visium_spatial(args.spatial, args.spatialo ,args.library_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spatial Drug Response Pipeline")

    parser.add_argument('--sc_tumor', type=str, default='data/sc-tumor/gse169246.h5ad')
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--preprocessed_sc', type=str, default='data/preprocessed/sc_tumor_preprocessed.h5ad',
                    help='Path to save preprocessed single-cell data')
    parser.add_argument('--spatial', type=str, default='data/spatial/visium-1142243F.h5ad',
                    help='Raw Visium spatial data file path')
    parser.add_argument('--spatialo', type=str, default='data/preprocessed/spatial_preprocessed.h5ad',
                    help='Path to save preprocessed spatial data')
    parser.add_argument('--library_id', type=str, default='1142243F', 
                       help='Library ID for the spatial data')
    args = parser.parse_args()




    # Preprocess data

    # Run main function with args object
    main(args)
