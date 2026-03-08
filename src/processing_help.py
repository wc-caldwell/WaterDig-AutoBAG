import os
from pygeohydro import EHydro
from concurrent.futures import ThreadPoolExecutor
from zipfile import ZipFile, BadZipFile

def extract_single_zip(zipf):
    """Extract a single zip file and remove it after extraction"""
    try:
        extract_dir = zipf.parent / zipf.stem
        extract_dir.mkdir(parents=True, exist_ok=True)
        with ZipFile(zipf, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        os.remove(zipf)
        return zipf.name, True
    except BadZipFile:
        os.remove(zipf)
        return zipf.name, False
    except Exception as e:
        return zipf.name, False

def retrieve_ehydro_data(data_dir, start_date, end_date, district_symbol, channel_area=None, surveyID = None, max_workers=8):
    """
    Function to retrieve eHydro survey data from the USACE database
    
    Parameters:
    -----------
    data_dir : Path
        Working directory for data storage
    start_date : str
        Start date in YYYY-MM-DD format
    end_date : str
        End date in YYYY-MM-DD format
    district_symbol : str
        USACE district code (e.g., 'CESWG')
    channel_area : str
        USACE NCF channel area named (e.g., 'CESWG_AN_01_BAY')
    surveyID : str
        USACE eHydro survey identifier (e.g., 'AN_01_BAY_20250523_CS')
    max_workers : int
        Number of parallel workers for zip extraction (default: 8)
    """
    where_clause = f"surveydatestart >= '{start_date}' AND surveydatestart <= '{end_date}' AND usacedistrictcode= '{district_symbol}'"
    data_dir = data_dir / district_symbol
    if surveyID:
        where_clause += f" AND surveyjobidpk= '{surveyID}'"

    if channel_area:
        where_clause += f" AND channelareaidfk= '{channel_area}'"
    
    data_dir.mkdir(parents=True, exist_ok=True)
    
    ehydro = EHydro(data_type="outlines", cache_dir=data_dir)
    topobathy = ehydro.bysql(where_clause)

    topobathy.to_parquet(data_dir / 'ehydro.parquet')
    print(f'eHydro survey data saved locally to {data_dir}')
    print(f'Survey metadata saved to {data_dir / "ehydro.parquet"}')

    # Parallel zip extraction for speed
    zip_files = list(data_dir.glob('*.ZIP')) + list(data_dir.glob('*.zip'))
    if zip_files:
        print(f'Extracting {len(zip_files)} zip files in parallel...')
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(extract_single_zip, zip_files))
        extracted = sum(1 for _, success in results if success)
        print(f'Successfully extracted {extracted}/{len(zip_files)} files')

    return topobathy